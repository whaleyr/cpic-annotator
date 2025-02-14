#!/usr/bin/env python3

__author__ = 'BinglanLi'

import sys
import json
import concurrent.futures
import pandas as pd
from glob import glob
from pathlib import Path
from utilities import extract_pcat_json, get_json_file_names
from concurrent.futures import ALL_COMPLETED
from timeit import default_timer as timer


if __name__ == "__main__":
    import argparse

    # describe the tool
    parser = argparse.ArgumentParser(description='Extract results from PharmCAT JSONs into a TSV-formatted file.')

    # list input arguments
    parser.add_argument("-i", "--input-dir", type=str, metavar='<dir>',
                        help="Directory path to PharmCAT JSONs.")
    parser.add_argument('-S', '--sample-file', type=str, metavar='<txt_file>',
                        help='(Optional) a file containing a list of sample IDs, one sample at a line.')
    parser.add_argument("-m", "--matcher-json-pattern", type=str, metavar='*.match.json',
                        default='*match.json',
                        help="(Optional) file name pattern of the Named Allele Matcher JSON files.")
    parser.add_argument("-p", "--phenotyper-json-pattern", type=str, metavar='*.phenotype.json',
                        default='*phenotype.json',
                        help="(Optional) file name pattern of the Phenotyper JSON files.")
    parser.add_argument("-a", "--allele-definition-files", type=str,
                        metavar='</path/to/*_translation.json>',
                        default=None,
                        help="(Optional) Path to PharmCAT allele definition JSON files. "
                             "Pattern matching strings should be included.")
    parser.add_argument("-gs", "--guideline-source", type=str, metavar='CPIC/DPWG',
                        default='CPIC',
                        help="(Optional) Guideline source to extract, default = CPIC.")
    parser.add_argument("-g", "--genes", type=str, metavar='gene1,gene2',
                        default=None,
                        help="(Optional) List of genes to be processed, separated by comma.")

    # output args
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument("-o", "--output-dir", type=str, metavar='<dir>',
                              help="(Optional) directory for outputs.  Defaults to the directory of the input.")
    output_group.add_argument("-bf", "--base-filename", type=str, metavar='<name>',
                              default='pharmcat_json2tsv',
                              help="(Optional) output prefix (without file extensions), "
                                   "by default \"pharmcat_json2tsv\".")

    # concurrency args
    concurrency_group = parser.add_argument_group('Concurrency arguments')
    concurrency_group.add_argument("-c", "--concurrent-mode", action="store_true",
                                   help="(Optional) use multiple processes - maximum number of processes spawned will "
                                        "default to to two less than the number of cpu cores.")
    concurrency_group.add_argument("-cp", "--max-concurrent-processes", type=int, metavar='<num processes>',
                                   default=None,
                                   help='(Optional) the maximum number of processes to use when concurrent mode ' +
                                        'is enabled.')

    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="(Optional) print more verbose messages")

    # parse arguments
    args = parser.parse_args()

    try:
        # parse args into internal variables
        m_input_dir: Path = Path(args.input_dir).resolve()
        if not m_input_dir.is_dir():
            print(f'-i/--input-dir is not a valid directory path: {m_input_dir}')
            sys.exit(1)
        print(f'Reading from the input directory {m_input_dir}')

        # define output file name
        m_output_basename: str = args.base_filename
        # define output file path
        if args.output_dir:
            m_output_dir: Path = Path(args.output_dir).resolve()
        else:
            m_output_dir: Path = Path.cwd()
            print(f'Output directory is not specified.\n'
                  f'\tSave to the current working directory:{m_output_dir}')
        # define the full output file path
        m_output_file: Path = m_output_dir / (m_output_basename + '.tsv')

        # get the list of matcher and phenotyper jsons
        m_matcher_jsons: list[str] = glob(str(m_input_dir.joinpath(args.matcher_json_pattern)))
        m_phenotyper_jsons: list[str] = glob(str(m_input_dir.joinpath(args.phenotyper_json_pattern)))
        # get sample list based on json file names
        matcher_json_samples: list[str] = ['.'.join(x.split('.')[:-2]) for x in m_matcher_jsons]
        phenotyper_json_samples: list[str] = ['.'.join(x.split('.')[:-2]) for x in m_phenotyper_jsons]
        # further narrow down sample list based on sample_file
        m_samples: list[str] = []
        if args.sample_file:
            # todo: replace this with preprocessor.util.read_sample_file
            print(f'Reading samples from file: {args.sample_file}')
            with open(args.sample_file, 'r') as f:
                # this assumes a small, non-memory-intensive sample file
                m_samples = [line.strip() for line in f.readlines()]
            json_samples: list[str] = list(set(matcher_json_samples) & set(phenotyper_json_samples) & set(m_samples))
        else:
            print('No sample file found. Reading all samples.')
            json_samples: list[str] = list(set(matcher_json_samples) & set(phenotyper_json_samples))
        # print out the total number of samples
        print(f'\tFound n = {len(json_samples)}')

        # specify genes to be processed
        m_genes: list[str] = []
        if args.guideline_source == 'CPIC':
            m_genes = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C9", "CYP2C19", "CYP3A5", "CYP4F2",
                       "DPYD", "G6PD", "IFNL3", "NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1",
                       "VKORC1", "CYP2D6", "HLA-A", "HLA-B", "MT-RNR1"]
        elif args.guideline_source == 'DPWG':
            m_genes = ["ABCG2", "CYP2B6", "CYP2C9", "CYP2C19", "CYP3A4", "CYP3A5", "DPYD", "NUDT15",
                       "SLCO1B1", "TPMT", "UGT1A1", "VKORC1", "CYP2D6"]
        # further narrow down the list of genes to be process if --genes is specified
        if args.genes is not None:
            m_selected_genes: list[str] = [x for x in args.genes.split(",") if x in m_genes]
            m_genes = m_selected_genes

        # read allele definition json files
        if args.allele_definition_files is None:
            script_dir: Path = Path(globals().get("__file__", "./_")).absolute().parent
            allele_definition_dir: Path = script_dir.parent.parent / ('main/resources/org/pharmgkb/pharmcat/definition'
                                                                      '/alleles')
            allele_definition_pattern: str = str(allele_definition_dir) + '/*_translation.json'
        else:
            m_allele_definition_jsons: Path = Path(args.allele_definition_files).resolve()
            if m_allele_definition_jsons.is_dir():
                allele_definition_pattern: str = str(m_allele_definition_jsons) + '/*_translation.json'
            else:
                allele_definition_pattern = str(m_allele_definition_jsons)
        print(f'Looking for the allele definition JSON files: {allele_definition_pattern}')
        allele_definition_jsons: list[str] = glob(allele_definition_pattern)
        # error out if no allele definition json is found
        if len(allele_definition_jsons) == 0:
            print('No allele definition JSON files found.\n')
            sys.exit(1)
        # read reference alleles at each allele-defining position
        allele_definition_references: dict[str, list[str]] = {}
        for json_file in allele_definition_jsons:
            print(f'\tReading {json_file}')
            # get the gene name
            gene: str = json_file.split('/')[-1].split('_')[0]
            # read file
            with open(json_file, 'r') as f:
                json_data: dict = json.load(f)
            # get reference genotype list
            positions: list[str] = [str(entry['position']) for entry in json_data['variants']]
            genotypes: list[str] = [entry['ref'] for entry in json_data['variants']]
            allele_definition_references[gene] = [x + ':' + y for x, y in zip(positions, genotypes)]

        # set up concurrent mode
        m_max_processes: int = 1
        if args.concurrent_mode:
            print('Concurrent mode enabled...')
        elif args.max_concurrent_processes is not None:
            print("-cp/--max_processes will be ignored (not running in multiprocess mode)")
        if args.max_concurrent_processes:
            m_max_processes = args.max_concurrent_processes

        start = timer()

        # extract json results for each sample
        summary_results: dict[str, list[str]] = {
            'sample': [], 'gene': [], 'phenotype': [], 'activity_score': [],
            'diplotype': [],
            'dpyd_ryr1_variants': [], 'dpyd_ryr1_variant_functions': [], 'dpyd_ryr1_variant_genotypes': [],
            'haplotype_1': [], 'haplotype_2': [],
            'haplotype_1_functions': [], 'haplotype_2_functions': [],
            'haplotype_1_variants': [], 'haplotype_2_variants': [],
            'missing_positions': [], 'uncallable_haplotypes': [], 
            'non_ref_genotypes': []
        }
        if args.concurrent_mode:
            with concurrent.futures.ProcessPoolExecutor(max_workers=m_max_processes) as e:
                futures = []
                print(f'Processing samples:')
                for one_sample in json_samples:
                    print(f'\t{one_sample}')

                    # get the matcher and phenotyper json names for the sample
                    matcher_file, phenotyper_file = get_json_file_names(one_sample, m_matcher_jsons, m_phenotyper_jsons)
                    futures.append(e.submit(extract_pcat_json, matcher_file, phenotyper_file, m_genes, one_sample,
                                            allele_definition_references, args.guideline_source))
                concurrent.futures.wait(futures, return_when=ALL_COMPLETED)
                for future in futures:
                    tmp_results = future.result()
                    # concatenate temporary dictionary with the summary dictionary
                    for key in summary_results:
                        summary_results[key].extend(tmp_results[key])
        else:
            print(f'Processing samples:')
            for one_sample in json_samples:
                print(f'\t{one_sample}')

                # get the matcher and phenotyper json names for the sample
                matcher_file, phenotyper_file = get_json_file_names(one_sample, m_matcher_jsons, m_phenotyper_jsons)

                # extract results for a sample
                tmp_results: dict[str, list[str]] = extract_pcat_json(
                    matcher_json=matcher_file,
                    phenotyper_json=phenotyper_file,
                    genes=m_genes,
                    sample_id=one_sample,
                    reference_genotypes=allele_definition_references,
                    guideline_source=args.guideline_source
                )

                # concatenate temporary dictionary with the summary dictionary
                for key in summary_results:
                    summary_results[key].extend(tmp_results[key])

        # convert the summary dictionary to a pandas data frame
        rez = pd.DataFrame.from_dict(summary_results)
        # remove null
        rez.loc[rez.phenotype == 'n/a', 'phenotype'] = ''
        # write to an output
        cols: list[str] = ['Sample', 'Gene', 'Phenotype', 'Activity_Score',
                           'Diplotype',
                           'DPYD_RYR1_Variants', 'DPYD_RYR1_Variant_Functions', 'DPYD_RYR1_Variant_Genotypes',
                           'Haplotype_1', 'Haplotype_2',
                           'Haplotype_1_Functions', 'Haplotype_2_Functions',
                           'Haplotype_1_Variants', 'Haplotype_2_Variants',
                           'Missing_Positions', 'Uncallable_Haplotypes', 'NonRefGenotypes']
        rez.to_csv(m_output_file.absolute(), mode='w', sep="\t", header=cols, index=False)

        end = timer()
        print("Done.")
        print("Preprocessed input VCF in %.2f seconds" % (end - start))

    except Exception as e:
        print(e)
        sys.exit(1)
