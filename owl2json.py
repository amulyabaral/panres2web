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
                # Save previous sequence
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Start new sequence
                # Extract ID from header (e.g., >pan_1 -> pan_1, >Pan_1_v1.0.1_identical -> Pan_1)
                header = line[1:].split()[0]  # Remove '>' and take first part
                # Remove version suffixes like _v1.0.1_identical
                if '_v' in header and '_identical' in header:
                    current_id = header.split('_v')[0]
                else:
                    current_id = header
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
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

    # Parse FASTA files
    print("\nParsing sequence files...")
    gene_sequences = parse_fasta(GENES_FASTA)
    protein_sequences = parse_fasta(PROTEINS_FASTA)

    # Data structure
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

    # First pass: collect all subjects and their properties
    print("Processing triples...")
    subject_props = defaultdict(lambda: {'types': [], 'properties': defaultdict(list)})

    # Track equivalent classes
    equivalences = {}  # Maps class -> canonical class

    for subj, pred, obj in graph:
        if not isinstance(subj, URIRef):
            continue  # Skip blank nodes

        subj_id = clean_uri(str(subj))
        pred_id = clean_uri(str(pred))

        if isinstance(obj, Literal):
            obj_value = str(obj)
            obj_data = {'value': obj_value, 'is_literal': True}
        elif isinstance(obj, URIRef):
            obj_id = clean_uri(str(obj))
            obj_data = {'value': obj_id, 'is_literal': False}
        else:
            continue  # Skip blank nodes

        # Track equivalentClass relationships
        if pred_id == 'owl:equivalentClass':
            # The object of the equivalentClass property is the canonical name
            # Pattern: Antibiotic_X equivalentClass Y -> Y is canonical
            canonical = obj_id
            alias = subj_id
            equivalences[alias] = canonical
            print(f"  Found equivalence: {alias} = {canonical}")

        # Store property
        subject_props[subj_id]['properties'][pred_id].append(obj_data)

        # Track types separately
        if pred_id == 'rdf:type':
            subject_props[subj_id]['types'].append(obj_data['value'])

    print(f"Found {len(subject_props)} subjects")
    print(f"Found {len(equivalences)} equivalence relationships")

    # Merge equivalent classes
    if equivalences:
        print("Merging equivalent classes...")
        for alias, canonical in equivalences.items():
            if alias in subject_props and canonical in subject_props:
                # Merge properties from alias into canonical
                for pred, values in subject_props[alias]['properties'].items():
                    if pred != 'owl:equivalentClass':  # Don't copy equivalence relationships
                        subject_props[canonical]['properties'][pred].extend(values)

                # Merge types
                for type_val in subject_props[alias]['types']:
                    if type_val not in subject_props[canonical]['types']:
                        subject_props[canonical]['types'].append(type_val)

                # Add alias as alternative label
                if 'rdfs:label' in subject_props[alias]['properties']:
                    alias_label = subject_props[alias]['properties']['rdfs:label'][0]['value']
                    subject_props[canonical]['properties']['alternative_labels'] = [
                        {'value': alias_label, 'is_literal': True}
                    ]

                # Remove the alias from subject_props
                del subject_props[alias]
                print(f"  Merged {alias} into {canonical}")

        # Update references to aliases to point to canonical (only for PanGenes)
        print("Updating references to equivalent classes in PanGenes...")
        for subj_id, props in subject_props.items():
            # Only update references for PanGenes, not OriginalGenes
            is_pangene = 'PanGene' in props['types'] or any(t in {'AntimicrobialResistanceGene', 'BiocideResistanceGene', 'MetalResistanceGene'} for t in props['types'])
            is_original_gene = 'OriginalGene' in props['types']

            if is_pangene and not is_original_gene:
                for pred, values in props['properties'].items():
                    for value in values:
                        if not value['is_literal'] and value['value'] in equivalences:
                            value['value'] = equivalences[value['value']]

    # Second pass: build categorized structure (OPTIMIZED)
    print("Categorizing subjects...")

    # Track which subjects are used as objects for categorization
    class_subjects = set()
    phenotype_subjects = set()
    mechanism_subjects = set()
    database_subjects = set()

    # Single pass through all subjects to collect relationships
    print(" -> Scanning for category relationships...")
    for subj_id, props in subject_props.items():
        if 'has_resistance_class' in props['properties']:
            for val in props['properties']['has_resistance_class']:
                if not val['is_literal']:
                    class_subjects.add(val['value'])

        if 'has_predicted_phenotype' in props['properties']:
            for val in props['properties']['has_predicted_phenotype']:
                if not val['is_literal']:
                    phenotype_subjects.add(val['value'])

        if 'has_mechanism_of_resistance' in props['properties']:
            for val in props['properties']['has_mechanism_of_resistance']:
                if not val['is_literal']:
                    mechanism_subjects.add(val['value'])

        if 'is_from_database' in props['properties']:
            for val in props['properties']['is_from_database']:
                if not val['is_literal']:
                    database_subjects.add(val['value'])

    print(f" -> Found {len(class_subjects)} classes, {len(phenotype_subjects)} phenotypes, "
          f"{len(mechanism_subjects)} mechanisms, {len(database_subjects)} databases")

    # Build final subject entries
    print(" -> Building subject entries...")
    total = len(subject_props)
    count = 0
    for subj_id, props in subject_props.items():
        count += 1
        if count % 1000 == 0:
            print(f"    Processed {count}/{total} subjects ({count*100//total}%)")

        # Get label
        label = subj_id
        if 'rdfs:label' in props['properties'] and props['properties']['rdfs:label']:
            label = props['properties']['rdfs:label'][0]['value']

        # Build subject entry
        subject_entry = {
            'id': subj_id,
            'label': label,
            'types': props['types'],
            'properties': {}
        }

        # Convert properties to simple format
        # For PanGenes, deduplicate property values (especially after equivalence mapping)
        is_pangene = 'PanGene' in props['types'] or any(t in {'AntimicrobialResistanceGene', 'BiocideResistanceGene', 'MetalResistanceGene'} for t in props['types'])

        for pred, values in props['properties'].items():
            if is_pangene:
                # Deduplicate values based on their content
                seen = set()
                deduped_values = []
                for val in values:
                    # Create a key that represents this value uniquely
                    if val['is_literal']:
                        key = ('literal', val['value'])
                    else:
                        key = ('ref', val['value'])

                    if key not in seen:
                        seen.add(key)
                        deduped_values.append(val)

                subject_entry['properties'][pred] = deduped_values
            else:
                # For non-PanGenes (like OriginalGenes), keep all values
                subject_entry['properties'][pred] = values

        # Add sequences if this is a gene or protein
        # Genes are in original format (various patterns from databases)
        # Proteins now use PAN* format (e.g., PAN1, PAN2) instead of Pan_* format
        if subj_id in gene_sequences:
            subject_entry['gene_sequence'] = gene_sequences[subj_id]
        if subj_id in protein_sequences:
            subject_entry['protein_sequence'] = protein_sequences[subj_id]

        data['subjects'][subj_id] = subject_entry

        # Categorize
        types = set(props['types'])

        # PanGene types
        pangene_types = {'PanGene', 'AntimicrobialResistanceGene', 'BiocideResistanceGene', 'MetalResistanceGene'}
        if types & pangene_types:
            data['categories']['PanGene'].append(subj_id)

        # OriginalGene
        if 'OriginalGene' in types:
            data['categories']['OriginalGene'].append(subj_id)

        # Role-based categories
        if subj_id in class_subjects:
            data['categories']['AntibioticClass'].append(subj_id)

        if subj_id in phenotype_subjects:
            data['categories']['Phenotype'].append(subj_id)

        if subj_id in mechanism_subjects:
            data['categories']['Mechanism'].append(subj_id)

        if subj_id in database_subjects:
            data['categories']['Database'].append(subj_id)

    data['metadata']['total_subjects'] = len(data['subjects'])

    # Print statistics
    print("\n=== Conversion Statistics ===")
    print(f"Total subjects: {data['metadata']['total_subjects']}")
    print(f"PanGenes: {len(data['categories']['PanGene'])}")
    print(f"Original Genes: {len(data['categories']['OriginalGene'])}")
    print(f"Antibiotic Classes: {len(data['categories']['AntibioticClass'])}")
    print(f"Phenotypes: {len(data['categories']['Phenotype'])}")
    print(f"Mechanisms: {len(data['categories']['Mechanism'])}")
    print(f"Databases: {len(data['categories']['Database'])}")

    # Write JSON
    print(f"\nWriting JSON to: {json_file}")
    with open(json_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    print(f"Done! JSON file created: {json_file}")
    print(f"File size: {len(json.dumps(data)) / 1024 / 1024:.2f} MB")

if __name__ == '__main__':
    # Download the latest OWL file from GitHub
    print("=== PanRes OWL to JSON Converter ===\n")

    if download_owl_file(GITHUB_RAW_URL, OWL_FILE):
        print()
        convert_owl_to_json(OWL_FILE, JSON_OUTPUT)
    else:
        print("\nFailed to download OWL file. Exiting.")
        sys.exit(1)
