#!/usr/bin/env python3
"""
Converts an OWL/RDF file to a simplified JSON format for use in static web applications.
Automatically fetches the latest OWL file from the PanResOntology GitHub repository.
"""

import json
import rdflib
from rdflib import URIRef, Literal
from collections import defaultdict
import sys
import urllib.request
import os

# Configuration
GITHUB_RAW_URL = 'https://raw.githubusercontent.com/rambiolab/PanResOntology/refs/heads/master/ontology/panres_v2.owl'
OWL_FILE = 'panres_v2.owl'
JSON_OUTPUT = 'panres2.json'
BASE_IRI = "http://genepi.dk/PanResOntology.owl#"
GENES_FASTA = 'panres2_genes.fa'
PROTEINS_FASTA = 'panres_final_protein.faa'

# Namespace prefixes
NAMESPACES = {
    BASE_IRI: '',
    "http://www.w3.org/2001/XMLSchema#": 'xsd:',
    "http://www.w3.org/1999/02/22-rdf-syntax-ns#": 'rdf:',
    "http://www.w3.org/2000/01/rdf-schema#": 'rdfs:',
    "http://www.w3.org/2002/07/owl#": 'owl:',
}


def download_owl_file(url, output_file):
    """
    Download the OWL file from GitHub repository.
    """
    print(f"Downloading OWL file from: {url}")
    try:
        with urllib.request.urlopen(url) as response:
            content = response.read()
            with open(output_file, 'wb') as f:
                f.write(content)
        print(f"Successfully downloaded to: {output_file}")
        print(f"File size: {len(content) / 1024:.2f} KB")
        return True
    except Exception as e:
        print(f"Error downloading OWL file: {e}")
        return False


def clean_uri(uri_str):
    """Convert full URI to shortened form."""
    for namespace, prefix in NAMESPACES.items():
        if uri_str.startswith(namespace):
            fragment = uri_str.split('#')[-1] if '#' in uri_str else uri_str.split('/')[-1]
            return prefix + fragment if prefix else fragment
    return uri_str


def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a dictionary mapping sequence IDs to sequences.
    """
    sequences = {}
    if not os.path.exists(fasta_file):
        print(f"Warning: FASTA file not found: {fasta_file}")
        return sequences

    print(f"Parsing FASTA file: {fasta_file}")
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)

                header = line[1:].split()[0]

                # Remove version suffixes like _v1.0.1_identical
                if '_v' in header and '_identical' in header:
                    current_id = header.split('_v')[0]
                else:
                    current_id = header

                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            sequences[current_id] = ''.join(current_seq)

    print(f"  Parsed {len(sequences)} sequences")
    return sequences


def convert_owl_to_json(owl_file, json_file):
    """
    Parse OWL file and create a JSON structure with:
    - subjects: dictionary of all subjects with their properties
    - categories: pre-categorized lists for quick access
    """
    print(f"Loading OWL file: {owl_file}")
    graph = rdflib.Graph()

    try:
        graph.parse(owl_file)
        print(f"Successfully parsed {len(graph)} triples")
    except Exception as e:
        print(f"Error parsing OWL file: {e}")
        sys.exit(1)

    print("\nParsing sequence files...")
    gene_sequences = parse_fasta(GENES_FASTA)
    protein_sequences = parse_fasta(PROTEINS_FASTA)

    data = {
        'subjects': {},
        'categories': {
            'PanGene': [],
            'OriginalGene': [],
            'AntibioticClass': [],
            'Phenotype': [],
            'Mechanism': [],
            'Database': []
        },
        'metadata': {
            'total_triples': len(graph),
            'total_subjects': 0
        }
    }

    print("Processing triples...")
    subject_props = defaultdict(lambda: {'types': [], 'properties': defaultdict(list)})

    # alias -> canonical/main class
    # Example:
    #   Antibiotic_Efflux owl:equivalentClass Efflux_Pump
    # means:
    #   keep Efflux_Pump, remove/canonicalize Antibiotic_Efflux
    equivalences = {}

    for subj, pred, obj in graph:
        if not isinstance(subj, URIRef):
            continue

        subj_id = clean_uri(str(subj))
        pred_id = clean_uri(str(pred))

        if isinstance(obj, Literal):
            obj_data = {'value': str(obj), 'is_literal': True}
        elif isinstance(obj, URIRef):
            obj_id = clean_uri(str(obj))
            obj_data = {'value': obj_id, 'is_literal': False}
        else:
            continue

        if pred_id == 'owl:equivalentClass' and not obj_data['is_literal']:
            alias = subj_id
            canonical = obj_data['value']
            equivalences[alias] = canonical
            print(f"  Found equivalence: {alias} -> {canonical}")

        subject_props[subj_id]['properties'][pred_id].append(obj_data)

        if pred_id == 'rdf:type':
            subject_props[subj_id]['types'].append(obj_data['value'])

    print(f"Found {len(subject_props)} subjects")
    print(f"Found {len(equivalences)} equivalence relationships")

    def canonical_id(entity_id):
        """
        Follow equivalence links until the final canonical/main ID.
        Protects against accidental cycles.
        """
        seen = set()
        current = entity_id

        while current in equivalences and current not in seen:
            seen.add(current)
            current = equivalences[current]

        return current

    # Collapse equivalence chains first, e.g. A -> B -> C becomes A -> C
    if equivalences:
        equivalences = {
            alias: canonical_id(canonical)
            for alias, canonical in equivalences.items()
        }

    def canonical_id(entity_id):
        """
        Follow collapsed equivalence links until the final canonical/main ID.
        """
        seen = set()
        current = entity_id

        while current in equivalences and current not in seen:
            seen.add(current)
            current = equivalences[current]

        return current

    # Merge alias subjects into canonical subjects, then remove aliases from output
    alias_subjects_to_skip = set()

    if equivalences:
        print("Merging equivalent classes into canonical/main classes...")

        for alias, canonical in list(equivalences.items()):
            if alias == canonical:
                continue

            alias_subjects_to_skip.add(alias)

            if alias in subject_props and canonical in subject_props:
                # Merge alias properties into canonical subject
                for pred, values in subject_props[alias]['properties'].items():
                    if pred == 'owl:equivalentClass':
                        continue
                    subject_props[canonical]['properties'][pred].extend(values)

                # Merge alias types into canonical subject
                for type_val in subject_props[alias]['types']:
                    if type_val not in subject_props[canonical]['types']:
                        subject_props[canonical]['types'].append(type_val)

                # Keep alias label as alternative label
                if 'rdfs:label' in subject_props[alias]['properties']:
                    alias_label = subject_props[alias]['properties']['rdfs:label'][0]['value']

                    existing_alt_labels = {
                        v['value']
                        for v in subject_props[canonical]['properties'].get('alternative_labels', [])
                        if v.get('is_literal')
                    }

                    if alias_label not in existing_alt_labels:
                        subject_props[canonical]['properties']['alternative_labels'].append({
                            'value': alias_label,
                            'is_literal': True
                        })

                print(f"  Merged {alias} into {canonical}")

            elif alias in subject_props and canonical not in subject_props:
                print(f"  Warning: canonical subject {canonical} not found for alias {alias}")

        # Rewrite ALL object references everywhere, not only PanGenes
        print("Updating references to equivalent classes everywhere...")
        for subj_id, props in subject_props.items():
            for pred, values in props['properties'].items():
                for value in values:
                    if not value['is_literal']:
                        value['value'] = canonical_id(value['value'])

    print("Categorizing subjects...")

    class_subjects = set()
    phenotype_subjects = set()
    mechanism_subjects = set()
    database_subjects = set()

    print(" -> Scanning for category relationships...")
    for subj_id, props in subject_props.items():
        if subj_id in alias_subjects_to_skip:
            continue

        if 'has_resistance_class' in props['properties']:
            for val in props['properties']['has_resistance_class']:
                if not val['is_literal']:
                    class_subjects.add(canonical_id(val['value']))

        if 'has_predicted_phenotype' in props['properties']:
            for val in props['properties']['has_predicted_phenotype']:
                if not val['is_literal']:
                    phenotype_subjects.add(canonical_id(val['value']))

        if 'has_mechanism_of_resistance' in props['properties']:
            for val in props['properties']['has_mechanism_of_resistance']:
                if not val['is_literal']:
                    mechanism_subjects.add(canonical_id(val['value']))

        if 'is_from_database' in props['properties']:
            for val in props['properties']['is_from_database']:
                if not val['is_literal']:
                    database_subjects.add(canonical_id(val['value']))

    print(f" -> Found {len(class_subjects)} classes, {len(phenotype_subjects)} phenotypes, "
          f"{len(mechanism_subjects)} mechanisms, {len(database_subjects)} databases")

    print(" -> Building subject entries...")
    total = len(subject_props)
    count = 0

    for subj_id, props in subject_props.items():
        count += 1

        if count % 1000 == 0:
            print(f"    Processed {count}/{total} subjects ({count * 100 // total}%)")

        # Do not write alias subjects as separate JSON subjects
        if subj_id in alias_subjects_to_skip:
            continue

        # Extra safety: if this subject canonicalizes to another one, skip it
        if canonical_id(subj_id) != subj_id:
            continue

        label = subj_id
        if 'rdfs:label' in props['properties'] and props['properties']['rdfs:label']:
            label = props['properties']['rdfs:label'][0]['value']

        subject_entry = {
            'id': subj_id,
            'label': label,
            'types': list(dict.fromkeys(props['types'])),
            'properties': {}
        }

        # Deduplicate all properties for all subjects
        for pred, values in props['properties'].items():
            seen = set()
            deduped_values = []

            for val in values:
                if not val['is_literal']:
                    val = {
                        'value': canonical_id(val['value']),
                        'is_literal': False
                    }

                if val['is_literal']:
                    key = ('literal', val['value'])
                else:
                    key = ('ref', val['value'])

                if key not in seen:
                    seen.add(key)
                    deduped_values.append(val)

            subject_entry['properties'][pred] = deduped_values

        if subj_id in gene_sequences:
            subject_entry['gene_sequence'] = gene_sequences[subj_id]

        if subj_id in protein_sequences:
            subject_entry['protein_sequence'] = protein_sequences[subj_id]

        data['subjects'][subj_id] = subject_entry

        types = set(props['types'])

        pangene_types = {
            'PanGene',
            'AntimicrobialResistanceGene',
            'BiocideResistanceGene',
            'MetalResistanceGene'
        }

        if types & pangene_types:
            data['categories']['PanGene'].append(subj_id)

        if 'OriginalGene' in types:
            data['categories']['OriginalGene'].append(subj_id)

        if subj_id in class_subjects:
            data['categories']['AntibioticClass'].append(subj_id)

        if subj_id in phenotype_subjects:
            data['categories']['Phenotype'].append(subj_id)

        if subj_id in mechanism_subjects:
            data['categories']['Mechanism'].append(subj_id)

        if subj_id in database_subjects:
            data['categories']['Database'].append(subj_id)

    # Remove duplicates and sort categories
    for cat in data['categories']:
        data['categories'][cat] = sorted(set(data['categories'][cat]))

    data['metadata']['total_subjects'] = len(data['subjects'])

    print("\n=== Conversion Statistics ===")
    print(f"Total subjects: {data['metadata']['total_subjects']}")
    print(f"PanGenes: {len(data['categories']['PanGene'])}")
    print(f"Original Genes: {len(data['categories']['OriginalGene'])}")
    print(f"Antibiotic Classes: {len(data['categories']['AntibioticClass'])}")
    print(f"Phenotypes: {len(data['categories']['Phenotype'])}")
    print(f"Mechanisms: {len(data['categories']['Mechanism'])}")
    print(f"Databases: {len(data['categories']['Database'])}")

    print("\nMechanisms kept in JSON:")
    for mech in data['categories']['Mechanism']:
        print(f"  {mech}")

    print(f"\nWriting JSON to: {json_file}")
    with open(json_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    print(f"Done! JSON file created: {json_file}")
    print(f"File size: {len(json.dumps(data)) / 1024 / 1024:.2f} MB")


if __name__ == '__main__':
    print("=== PanRes OWL to JSON Converter ===\n")

    if download_owl_file(GITHUB_RAW_URL, OWL_FILE):
        print()
        convert_owl_to_json(OWL_FILE, JSON_OUTPUT)
    else:
        print("\nFailed to download OWL file. Exiting.")
        sys.exit(1)
