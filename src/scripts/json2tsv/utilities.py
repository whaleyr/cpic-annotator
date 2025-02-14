#!/usr/bin/env python3

__author__ = 'BinglanLi'

import json
import pandas as pd


def get_json_file_names(sample_id: str, matcher_json_list: list[str], phenotyper_json_list: list[str]) \
        -> tuple[str, str]:
    # get the matcher and phenotyper json names for the sample
    matcher_pattern: str = f"{sample_id}.match.json"
    matcher_match_list: list[str] = [x for x in matcher_json_list if matcher_pattern in x]
    phenotyper_pattern: str = f"{sample_id}.phenotype.json"
    phenotyper_match_list: list[str] = [x for x in phenotyper_json_list if phenotyper_pattern in x]

    matcher_file: str = ''
    phenotyper_file: str = ''
    if len(matcher_match_list) == 1 and len(phenotyper_match_list) == 1:
        matcher_file = matcher_match_list[0]
        phenotyper_file = phenotyper_match_list[0]
    else:
        print(f'WARNING: multiple matches found for {sample_id}')

    return matcher_file, phenotyper_file


def get_names_and_genotypes(json_entry: dict, reference_genotypes: list[str]) -> tuple[list[str], list[str]]:

    # set up empty return variables
    haplotype_names: list[str] = []
    alternative_genotypes: list[str] = []

    if json_entry:
        # identify whether the json entry has component haplotypes
        has_component_haplotypes: bool = True if 'componentHaplotypes' in json_entry else False

        # get the list of alternative genotypes
        tmp_genotypes: set[str] = {y for x in json_entry['sequences'] for y in x.split(';')
                                   if y != '' and y not in reference_genotypes}

        # append variant/haplotype names and alternative genotypes to return variables
        if has_component_haplotypes:
            # get names
            for component_entry in json_entry['componentHaplotypes']:
                haplotype_names.append(component_entry['name'])
            # append alternative genotypes
            for i in range(len(haplotype_names)):
                if i == 0:
                    alternative_genotypes.append(','.join(tmp_genotypes))
                else:
                    alternative_genotypes.append('')
        else:
            # get names
            haplotype_names.append(json_entry['name'])
            # append alternative genotypes
            alternative_genotypes.append(','.join(tmp_genotypes))
    else:
        haplotype_names.append('')
        alternative_genotypes.append('')

    return haplotype_names, alternative_genotypes

def get_non_ref_genotypes(variants):
    non_ref_genotypes = []
    for variant in variants:
        wildtypeAllele = variant['wildtypeAllele']
        wildtypeGenotype = f"{wildtypeAllele}/{wildtypeAllele}"
        if variant['call'] != wildtypeGenotype:
            non_ref_genotypes.append(variant)
    
    non_ref_genotypes = sorted(non_ref_genotypes, key=lambda x: x['position'])
    non_ref_genotypes = [f"{x['position']}:{x['call']}" for x in non_ref_genotypes]

    return (';'.join(non_ref_genotypes))


def extract_pcat_json(matcher_json: str, phenotyper_json: str, genes: list[str], sample_id: str,
                      reference_genotypes: dict[str, list[str]], guideline_source: str = 'CPIC') \
        -> dict[str, list[str]]:
    # read json file
    with open(matcher_json, 'r') as mf:
        matcher_data: dict = json.load(mf)['results']
    with open(phenotyper_json, 'r') as pf:
        phenotyper_data: dict = json.load(pf)['geneReports'][guideline_source]

    # special toxicity genes
    special_toxicity_genes: list[str] = ['DPYD', 'RYR1']

    # get the gene list from the named allele matcher
    # to locate the index of a gene when extracting named allele matcher data
    m_gene_list: list[str] = [x['gene'] for x in matcher_data]

    # set up the summary dictionary
    pcat_summary: dict[str, list[str]] = {
        'sample': [], 'gene': [], 'phenotype': [], 'activity_score': [],
        'diplotype': [],
        'dpyd_ryr1_variants': [], 'dpyd_ryr1_variant_functions': [], 'dpyd_ryr1_variant_genotypes': [],
        'haplotype_1': [], 'haplotype_2': [],
        'haplotype_1_functions': [], 'haplotype_2_functions': [],
        'haplotype_1_variants': [], 'haplotype_2_variants': [],
        'missing_positions': [], 'uncallable_haplotypes': [], 
        'non_ref_genotypes': []
    }

    for i_gene in genes:
        # get the index of the gene in the named allele matcher json
        if i_gene in m_gene_list:
            idx_ma_gene = m_gene_list.index(i_gene)
        else:
            continue

        # get uncallable haplotypes
        uncallable_haplotypes: str = ';'.join(matcher_data[idx_ma_gene]['uncallableHaplotypes'])  # named allele matcher
        # get missing positions
        missing_positions: str = ';'.join(
            [str(x['position']) for x in matcher_data[idx_ma_gene]['matchData']['missingPositions']]
        )  # named allele matcher

        # get the number of diplotype calls for a gene
        # for special toxicity genes, like DPYD and RYR1, the number is always one
        # for the rest of the genes, the number is determined by the number of diplotypes in the matcher json
        n_diploids: int = 1 if i_gene in special_toxicity_genes else len(matcher_data[idx_ma_gene]['diplotypes'])

        # fill up the summary dictionary based on the number of PharmCAT calls
        if n_diploids == 0:
            # append empty entries to the summary dictionary
            pcat_summary['sample'].append(sample_id)
            pcat_summary['gene'].append(i_gene)
            pcat_summary['phenotype'].append('')
            pcat_summary['activity_score'].append('')
            pcat_summary['diplotype'].append('')
            pcat_summary['dpyd_ryr1_variants'].append('')
            pcat_summary['dpyd_ryr1_variant_functions'].append('')
            pcat_summary['dpyd_ryr1_variant_genotypes'].append('')
            pcat_summary['haplotype_1'].append('')
            pcat_summary['haplotype_2'].append('')
            pcat_summary['haplotype_1_functions'].append('')
            pcat_summary['haplotype_2_functions'].append('')
            pcat_summary['haplotype_1_variants'].append('')
            pcat_summary['haplotype_2_variants'].append('')
            pcat_summary['missing_positions'].append(missing_positions)
            pcat_summary['uncallable_haplotypes'].append(uncallable_haplotypes)
            pcat_summary['non_ref_genotypes'].append(get_non_ref_genotypes(phenotyper_data[i_gene]['variants']))
        elif i_gene in special_toxicity_genes:
            # get the phenotyper entries
            # in phenotyper json, special toxicity genes have only one entry of recommendationDiplotypes
            p_entry: dict = phenotyper_data[i_gene]['recommendationDiplotypes'][0]
            # append the extracted data into the summary dictionary
            pcat_summary['sample'].append(sample_id)
            pcat_summary['gene'].append(i_gene)
            pcat_summary['phenotype'].append(p_entry['phenotypes'][0] if len(p_entry['phenotypes']) else '')
            pcat_summary['activity_score'].append(
                str(p_entry['activityScore']) if p_entry['activityScore'] is not None else '')

            # identify whether the sample is
            # (1) effectively phased with component haplotypes,
            # (2) effectively phased without component haplotypes, or
            # (3) not effectively phased
            is_effectively_phased = matcher_data[idx_ma_gene]['matchData']['effectivelyPhased']

            # specify different json fields to extract for effectively phased and not effectively phased samples
            if is_effectively_phased:
                # append effectively phased diplotype
                pcat_summary['diplotype'].append(
                    phenotyper_data[i_gene]['sourceDiplotypes'][0]['label'])
                # specify json fields for effectively phased samples
                p_toxicity_variant_field: str = 'matcherComponentHaplotypes'
                # use the 'diplotypes', instead of the 'haplotypes' field to extract homozygous variants
                ma_toxicity_variant_field: str = 'diplotypes'
            else:
                # do not append diplotype for not effectively phased samples
                pcat_summary['diplotype'].append('')
                # specify json fields for not effectively phased samples
                p_toxicity_variant_field: str = 'sourceDiplotypes'
                ma_toxicity_variant_field: str = 'haplotypeMatches'

            # set up variables to collect individual toxicity variants
            p_toxicity_variants: list[str] = []
            p_toxicity_functions: list[str] = []
            ma_toxicity_variants: list[str] = []
            ma_toxicity_genotypes: list[str] = []

            # get component haplotypes and functions from the phenotyper json
            for entry in phenotyper_data[i_gene][p_toxicity_variant_field]:
                p_toxicity_variants.append(entry['allele1']['name'])
                p_toxicity_functions.append(entry['allele1']['function'])

            # get component haplotypes and genotypes from matcher json
            # for effectively phased diplotypes
            if is_effectively_phased:
                # get the toxicity variant name
                for json_field in ['haplotype1', 'haplotype2']:
                    # get the json field of a specific haplotype in a diplotype
                    ma_entry: dict = matcher_data[idx_ma_gene][ma_toxicity_variant_field][0][json_field]

                    # get variant names and alternative genotypes
                    tmp_names, tmp_genotypes = get_names_and_genotypes(ma_entry, reference_genotypes[i_gene])
                    ma_toxicity_variants.extend(tmp_names)
                    ma_toxicity_genotypes.extend(tmp_genotypes)
            else:  # for not effectively phased diplotypes
                for entry in matcher_data[idx_ma_gene][ma_toxicity_variant_field]:
                    tmp_names, tmp_genotypes = get_names_and_genotypes(entry, reference_genotypes[i_gene])
                    ma_toxicity_variants.extend(tmp_names)
                    ma_toxicity_genotypes.extend(tmp_genotypes)

            # concatenate individual toxicity variants into one entry
            df_p_toxicity_functions: pd.DataFrame = pd.DataFrame({
                "variants": p_toxicity_variants,
                "functions": p_toxicity_functions
            })
            df_ma_toxicity_genotypes: pd.DataFrame = pd.DataFrame({
                "variants": ma_toxicity_variants,
                "genotypes": ma_toxicity_genotypes
            })
            df_toxicity: pd.DataFrame = pd.merge(df_ma_toxicity_genotypes, df_p_toxicity_functions,
                                                 how='left', on="variants")

            # append the extracted data into the summary dictionary
            pcat_summary['dpyd_ryr1_variants'].append(';'.join(df_toxicity.variants))
            pcat_summary['dpyd_ryr1_variant_functions'].append(';'.join(df_toxicity.functions))
            pcat_summary['dpyd_ryr1_variant_genotypes'].append(
                ','.join(df_toxicity.genotypes[df_toxicity.genotypes != ''])
            )
            pcat_summary['haplotype_1'].append('')
            pcat_summary['haplotype_2'].append('')
            pcat_summary['haplotype_1_functions'].append('')
            pcat_summary['haplotype_2_functions'].append('')
            pcat_summary['haplotype_1_variants'].append('')
            pcat_summary['haplotype_2_variants'].append('')
            pcat_summary['missing_positions'].append(missing_positions)
            pcat_summary['uncallable_haplotypes'].append(uncallable_haplotypes)
            pcat_summary['non_ref_genotypes'].append(get_non_ref_genotypes(phenotyper_data[i_gene]['variants']))

        else:  # for genes other than DPYD and RYR1
            for idx_diplotype in range(n_diploids):
                # get the matcher and phenotyper entries
                # matcher and phenotyper jsons have sorted and listed diplotypes in the same order
                # todo: remove the if condition when duplicate cyp2d6 diplotypes are fixed in the research mode
                if len(phenotyper_data[i_gene]['recommendationDiplotypes']) - 1 >= idx_diplotype:
                    p_entry: dict = phenotyper_data[i_gene]['recommendationDiplotypes'][idx_diplotype]
                else:
                    continue
                ma_entry: dict = matcher_data[idx_ma_gene]['diplotypes'][idx_diplotype]

                # append the extracted phenotyper data into the summary dictionary
                pcat_summary['sample'].append(sample_id)
                pcat_summary['gene'].append(i_gene)
                pcat_summary['phenotype'].append(p_entry['phenotypes'][0] if len(p_entry['phenotypes']) else '')
                pcat_summary['activity_score'].append(
                    str(p_entry['activityScore']) if p_entry['activityScore'] is not None else '')
                pcat_summary['diplotype'].append(p_entry['label'])
                pcat_summary['dpyd_ryr1_variants'].append('')
                pcat_summary['dpyd_ryr1_variant_functions'].append('')
                pcat_summary['dpyd_ryr1_variant_genotypes'].append('')
                pcat_summary['haplotype_1'].append(p_entry['allele1']['name'])
                pcat_summary['haplotype_1_functions'].append(p_entry['allele1']['function'])
                # account for empty allele2 due to, for example, haploids in G6PD
                if p_entry['allele2']:
                    pcat_summary['haplotype_2'].append(p_entry['allele2']['name'])
                    pcat_summary['haplotype_2_functions'].append(p_entry['allele2']['function'])
                else:
                    pcat_summary['haplotype_2'].append('')
                    pcat_summary['haplotype_2_functions'].append('')
                pcat_summary['missing_positions'].append(missing_positions)
                pcat_summary['uncallable_haplotypes'].append(uncallable_haplotypes)

                # append the extracted matcher data into the summary dictionary
                for json_field in ['haplotype1', 'haplotype2']:
                    # get the alternative genotypes for a haplotype
                    tmp_names, tmp_genotypes = get_names_and_genotypes(
                        ma_entry[json_field], reference_genotypes[i_gene])
                    # append the alternative genotypes to the variant list
                    if json_field == 'haplotype1':
                        pcat_summary['haplotype_1_variants'].append(','.join(tmp_genotypes))
                    if json_field == 'haplotype2':
                        pcat_summary['haplotype_2_variants'].append(','.join(tmp_genotypes))
        
                pcat_summary['non_ref_genotypes'].append(get_non_ref_genotypes(phenotyper_data[i_gene]['variants']))
                
    return pcat_summary
