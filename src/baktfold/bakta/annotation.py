
import re

from typing import Sequence
from loguru import logger

import baktfold.bakta.config as cfg
import baktfold.bakta.constants as bc
import baktfold.io.insdc as insdc


RE_MULTIWHITESPACE = re.compile(r'\s{2,}')

RE_PROTEIN_CONTIG = re.compile(r'\s+contig\s*', flags=re.IGNORECASE)
RE_PROTEIN_HOMOLOG = re.compile(r'\shomolog(?: (\d+))?', flags=re.IGNORECASE)
RE_PROTEIN_PUTATIVE = re.compile(r'(potential|possible|probable|predicted)', flags=re.IGNORECASE)
RE_PROTEIN_NODE = re.compile(r'NODE_', flags=re.IGNORECASE)
RE_PROTEIN_POTENTIAL_CONTIG_NAME = re.compile(r'(genome|shotgun)', flags=re.IGNORECASE)
RE_PROTEIN_DOMAIN_CONTAINING = re.compile(r'domain-containing protein', flags=re.IGNORECASE)
RE_PROTEIN_REMNANT = re.compile(r'Remnant of ', re.IGNORECASE)
RE_PROTEIN_TMRNA = re.compile(r'TmRNA', flags=re.IGNORECASE)
RE_PROTEIN_NO_LETTERS = re.compile(r'[^A-Za-z]')
RE_PROTEIN_SUSPECT_CHARS_DISCARD = re.compile(r'[.#]')
RE_PROTEIN_SUSPECT_CHARS_REPLACE = re.compile(r'[@=?%]')
RE_PROTEIN_SUSPECT_CHARS_BEGINNING = r"_\-+.:,;/\\'"
RE_PROTEIN_PERIOD_SEPARATOR = re.compile(r'([a-zA-Z0-9]+)\.([a-zA-Z0-9]+)')
RE_PROTEIN_WRONG_PRIMES = re.compile(r'[\u2032\u0060\u00B4]')  # prime (′), grave accent (`), acute accent (´)
RE_PROTEIN_WEIGHT = re.compile(r' [0-9]+(?:\.[0-9]+)? k?da ', flags=re.IGNORECASE)
RE_PROTEIN_SYMBOL = re.compile(r'[A-Z][a-z]{2}[A-Z][0-9]?')
RE_DOMAIN_OF_UNKNOWN_FUNCTION = re.compile(r'(DUF\d{3,4})', flags=re.IGNORECASE)
RE_UNCHARACTERIZED_PROTEIN_FAMILY = re.compile(r'(UPF\d{3,4})', flags=re.IGNORECASE)
RE_GENE_CAPITALIZED = re.compile(r'^[A-Z].+', flags=re.DOTALL)
RE_GENE_SUSPECT_CHARS = re.compile(r'[\?]', flags=re.DOTALL)
RE_GENE_SYMBOL = re.compile(r'[a-z]{3}[A-Z][0-9]?')


def combine_annotation(feature: dict, fast: bool):
    """
    Combines annotation information from different sources into a single feature.

    Args:
      feature (dict): The feature to combine annotation for.
      fast (bool): If True, skips AFDB
    Returns:
      None

    Examples:
      >>> combine_annotation(feature)
    """


    # ups = feature.get('ups', None)
    # ips = feature.get('ips', None)
    # psc = feature.get('psc', None)
    # pscc = feature.get('pscc', None)
    pstc = feature.get('pstc', None)
    # expert_hits = feature.get('expert', [])

    # gene = None
    # genes = set()
    # product = None

    product = feature.get('product', None)
    db_xrefs = feature.get('db_xrefs', [])

    if(pstc):

        # Always normalize pstc to a list
        if isinstance(pstc, dict):
            pstc = [pstc]
        elif isinstance(pstc, str):
            pstc = [pstc]
        
        # afdb
        afdb_entry = None if fast else next(
            (p for p in pstc if isinstance(p, dict) and p.get('source') == 'afdb'),
            None
        )
        # swissprot
        swissprot_entry = next((p for p in pstc if isinstance(p, dict) and p.get('source') == 'swissprot'), None)
        # pdb
        pdb_entry = next((p for p in pstc if isinstance(p, dict) and p.get('source') == 'pdb'), None)
        # cath
        cath_entry = next((p for p in pstc if isinstance(p, dict) and p.get('source') == 'cath'), None)
        # custom
        custom_entry = next((p for p in pstc if isinstance(p, dict) and p.get('source') == 'custom_db'), None)

        ####
        # hierarchy
        # if it exists, custom is at the top
        # custom
        # if not
        # 1. SwissProt
        # 2. AFDB
        # 3. PDB
        # 4. CATH
        ####

        if custom_entry:
            pstc_product = custom_entry['description'] 
        elif swissprot_entry:
            pstc_product = swissprot_entry['description']
        elif afdb_entry:
            pstc_product = afdb_entry['description'] 
        elif pdb_entry:
            pstc_product = pdb_entry['description'] 
        elif cath_entry:
            pstc_product = cath_entry['description'] 
        else:
            pstc_product = None

        if(pstc_product):
            product = pstc_product

        # Collect all db_xref IDs
        for entry in pstc:
            if isinstance(entry, dict):
                src = entry.get('source', '').lower()
                eid = entry.get('id')
                if eid:
                    if src == 'afdb':
                        if not fast:
                            db_xrefs.append(f"afdb_v6:afdbclusters_{eid}")
                    elif src == 'swissprot':
                        db_xrefs.append(f"afdb_v6:swissprot_{eid}")
                    elif src == 'pdb':
                        db_xrefs.append(f"pdb:pdb_{eid}")
                    elif src == 'cath':
                        db_xrefs.append(f"cath:cath_{eid}")
                    elif src == 'custom_db':
                        db_xrefs.append(f"custom:custom_{eid}")
                    else:
                        db_xrefs.append(eid)
            elif isinstance(entry, str):
                # Preserve any existing string cross-references
                db_xrefs.append(entry)

        # mark as baktfold
        mark_as_baktfold(feature)




    # if(len(expert_hits) > 0):
    #     top_expert_hit = sorted(expert_hits,key=lambda k: (k['rank'], k.get('score', 0), calc_annotation_score(k)), reverse=True)[0]
    #     expert_genes = top_expert_hit.get('gene', None)
    #     if(expert_genes):
    #         expert_genes = expert_genes.replace('/', ',').split(',')
    #         genes.update(expert_genes)
    #         gene = expert_genes[0]
    #     product = top_expert_hit.get('product', None)
    #     for hit in expert_hits:
    #         db_xrefs.update(hit.get('db_xrefs', []))

    if product and "hypothetical protein" not in product.lower():
        product = revise_cds_product(product)
        if(product):
            if(cfg.compliant):
                product = insdc.revise_product_insdc(product)
            feature['product'] = product

            unmark_as_hypothetical(feature)
        
            # protein_gene_symbol = extract_protein_gene_symbol(product)
            # if(protein_gene_symbol):
            #     genes.add(protein_gene_symbol)
            # revised_genes = revise_cds_gene_symbols(genes)
            # revised_gene = None
            # if gene is not None:
            #     revised_gene = revise_cds_gene_symbols([gene])  # special treatment for selected gene symbol
            #     revised_gene = revised_gene[0] if len(revised_gene) > 0 else None
            # if(revised_gene is None  and  len(revised_genes) >= 1):  # select first from gene symbol list if no symbol was selected before
            #     revised_gene = revised_genes[0]

            # feature['gene'] = revised_gene
            # feature['genes'] = sorted(revised_genes)
        else:
            mark_as_hypothetical(feature)
    else:
        mark_as_hypothetical(feature)

    feature['db_xrefs'] = sorted(list(db_xrefs))




def calc_annotation_score(orf:dict) -> int:
    """
    Calculates the annotation score for a given ORF.

    Args:
      orf (dict): The ORF to calculate the annotation score for.

    Returns:
      int: The annotation score for the given ORF.

    Examples:
      >>> calc_annotation_score(orf)
    """
    score = 0
    if(orf.get('gene', None)):
        score += 1
    if(orf.get('product', None)):
        score += 1
    return score


def extract_protein_gene_symbol(product: str) -> str:
    """
    Extracts a valid gene symbol from a protein name.

    Args:
      product (str): The protein name to extract a gene symbol from.

    Returns:
      str: The extracted gene symbol.

    Examples:
      >>> extract_protein_gene_symbol(product)
    """
    gene_symbols = []
    for part in product.split(' '):  # try to extract valid gene symbols
        m = RE_GENE_SYMBOL.fullmatch(part)
        if(m):
            symbol = m[0]
            logger.info('fix gene: extract symbol from protein name. symbol=%s', symbol)
            gene_symbols.append(symbol)
        else:
            m = RE_PROTEIN_SYMBOL.fullmatch(part)  # extract protein names
            if(m):
                symbol = m[0]
                symbol = symbol[0].lower() + symbol[1:]
                logger.info('fix gene: extract symbol from protein name. symbol=%s', symbol)
                gene_symbols.append(symbol)
    if(len(gene_symbols) == 0):  # None found
        return None
    elif(len(gene_symbols) == 1):  # found 1
        return gene_symbols[0]
    else:  # found more than one, take the 2nd as the 1st often describes a broader gene family like "xyz family trancsriptional regulator ..."
        return gene_symbols[1]


def revise_cds_gene_symbols(raw_genes: Sequence[str]):
    """
    Revises a list of gene symbols to ensure they are valid.

    Args:
      raw_genes (Sequence[str]): The list of gene symbols to revise.

    Returns:
      list: The revised list of gene symbols.

    Examples:
      >>> revise_cds_gene_symbols(raw_genes)
    """
    revised_genes = set()
    for gene in raw_genes:
        old_gene = gene
        if(RE_GENE_SUSPECT_CHARS.search(gene)):  # check for suspect characters -> remove gene symbol
            logger.info('fix gene: remove gene symbol containing suspect chars. old=%s', old_gene)
            continue

        old_gene = gene
        gene = gene.replace('gene', '')
        if(gene != old_gene):  # remove gene literal
            logger.info('fix gene: remove gene literal. new=%s, old=%s', gene, old_gene)

        old_gene = gene
        if(gene[-1] == '-'):  # remove orphan hyphen
            gene = gene[:-1]
            logger.info('fix gene: remove orphan hypen. new=%s, old=%s', gene, old_gene)
        
        old_gene = gene
        gene = RE_MULTIWHITESPACE.sub(' ', gene).strip()  # revise whitespaces
        if(gene != old_gene):
            logger.info('fix gene: revise whitespaces. new=%s, old=%s', gene, old_gene)

        old_gene = gene
        if(RE_GENE_CAPITALIZED.fullmatch(gene)):
            gene = gene[0].lower() + gene[1:]
            logger.info('fix gene: lowercase first char. new=%s, old=%s', gene, old_gene)

        if(len(gene) >= 3):
            if(len(gene) <= 12):
                revised_genes.add(gene)
            else:
                old_gene = gene
                gene = extract_protein_gene_symbol(gene)
                if(gene):
                    revised_genes.add(gene)
    return list(revised_genes)


def revise_cds_product(product: str):
    """Revise product name for INSDC compliant submissions"""

    # from gb 
    # grep "Uncharacterized protein" AFDBClusters.tsv | wc -l
    #     805448

    if "Uncharacterized protein" in product:
        old_product = product
        product = "hypothetical protein"
        if product != old_product:
            logger.info(f'fix product: renamed uncharacterized protein as hypothetical. new={product}, old={old_product}')

    # from bakta

    old_product = product
    product = RE_PROTEIN_WEIGHT.sub(' ', product)  # remove protein weight in (k)Da
    if(product != old_product):
        logger.info('fix product: remove protein weight in (k)Da. new=%s, old=%s', product, old_product)

    old_product = product
    product = re.sub(RE_PROTEIN_PERIOD_SEPARATOR, r'\1-\2', product)  # replace separator periods
    if(product != old_product):
        logger.info('fix product: replace separator periods. new=%s, old=%s', product, old_product)
    
    old_product = product
    if(product[0] in RE_PROTEIN_SUSPECT_CHARS_BEGINNING):  # remove suspect first character
        product = product[1:]
        logger.info('fix product: replace invalid first character. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_SUSPECT_CHARS_DISCARD.sub('', product)  # remove suspect characters
    if(product != old_product):
        logger.info('fix product: replace invalid characters. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_SUSPECT_CHARS_REPLACE.sub(' ', product)  # replace suspect characters by single whitespace
    if(product != old_product):
        logger.info('fix product: replace invalid characters. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_WRONG_PRIMES.sub('\u0027', product)  # replace wrong prime characters with single quote (U+0027) (') according to https://www.ncbi.nlm.nih.gov/genome/doc/internatprot_nomenguide/
    if(product != old_product):
        logger.info('fix product: replace wrong prime characters. new=%s, old=%s', product, old_product)

    old_product = product
    product = product.replace('FOG:', '')  # remove FOG ids
    if(product != old_product):
        logger.info('fix product: replace FOG ids. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_REMNANT.sub('', product)  # remove 'Remnant of's
    if(product != old_product):
        logger.info('fix product: replace remnant ofs. new=%s, old=%s', product, old_product)

    old_product = product
    dufs = []  # replace DUF-containing products
    for m in RE_DOMAIN_OF_UNKNOWN_FUNCTION.finditer(product):
        dufs.append(m.group(1).upper())
    if(len(dufs) >= 1):
        product = f"{' '.join(dufs)} domain{'s' if len(dufs) > 1 else ''}-containing protein"
        if(product != old_product):
            logger.info('fix product: revise DUF. new=%s, old=%s', product, old_product)
    
    old_product = product
    if('conserved' in product.lower()):  # replace conserved UPF proteins
        upfs = []
        for m in RE_UNCHARACTERIZED_PROTEIN_FAMILY.finditer(product):
            upfs.append(m.group(1).upper())
        if(len(upfs) >= 1):
            product = f"{' '.join(upfs)} protein"
            if(product != old_product):
                logger.info('fix product: revise UPF. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_HOMOLOG.sub('-like protein', product)  # replace Homologs
    if(product != old_product):
        if(product.count('protein') == 2):
            product = product.replace('protein', '', 1)  # remove former protein term if existing
        logger.info('fix product: replace Homolog. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_MULTIWHITESPACE.sub(' ', product).strip()  # revise whitespaces
    if(product != old_product):
        logger.info('fix product: revise whitespaces. new=%s, old=%s', product, old_product)

    old_product = product
    product = RE_PROTEIN_PUTATIVE.sub('putative', product)  # replace putative synonyms)
    if(product != old_product):
        logger.info('fix product: replace putative synonyms. new=%s, old=%s', product, old_product)

    old_product = product
    if(RE_PROTEIN_DOMAIN_CONTAINING.search(product)):  # replace domain name underscores in domain names
        product = product.replace('_', '-')
        if(product != old_product):
            logger.info('fix product: replace domain name underscores. new=%s, old=%s', product, old_product)
    
    old_product = product
    if(RE_PROTEIN_TMRNA.fullmatch(product)):
        product = ''
        logger.info('fix product: discard pure tmRNA product descriptions. new=%s, old=%s', product, old_product)

    old_product = product
    if(
        RE_PROTEIN_CONTIG.search(product) or  # protein containing 'sequence'
        RE_PROTEIN_NODE.search(product) or  # potential contig name (SPAdes)
        RE_PROTEIN_POTENTIAL_CONTIG_NAME.search(product) or  # potential contig name (SPAdes)
        RE_PROTEIN_NO_LETTERS.fullmatch(product)  # no letters -> set to Hypothetical
        ):  # remove suspect products and mark as hypothetical
        product = None
        logger.info('remove product: mark proteins with suspect products as hypothetical. old=%s', old_product)

    return product


def mark_as_hypothetical(feature: dict):
    """
    Marks a feature as hypothetical.

    Args:
      feature (dict): The feature to mark as hypothetical.

    Returns:
      None

    Examples:
      >>> mark_as_hypothetical(feature)
    """
    # no need to actually print this I think
    # logger.info(
    #     f'marked as hypothetical: seq={feature['sequence']}, start={feature['start']}, stop={feature['stop']}, strand={feature['strand']}'
    # )
    feature['hypothetical'] = True
    feature['gene'] = None
    feature['genes'] = []
    feature['product'] = bc.HYPOTHETICAL_PROTEIN

def unmark_as_hypothetical(feature: dict):
    """
    Removes the hypothetical key from the given feature dictionary.

    Args:
      feature (dict): The feature dictionary to remove the hypothetical key from.

    Returns:
      None

    Examples:
      >>> feature = {'sequence': 'ATG', 'start': 1, 'stop': 3, 'strand': '+'}
      >>> unmark_as_hypothetical(feature)
      >>> feature
      {'sequence': 'ATG', 'start': 1, 'stop': 3, 'strand': '+'}
    """
    # logger.info(
    #     f'unmarked as hypothetical: seq={feature['sequence']}, start={feature['start']}, stop={feature['stop']}, strand={feature['strand']}'
    # )
    feature.pop('hypothetical', None)  # remove completely

def mark_as_baktfold(feature: dict):
    """
    Adds the baktfold key to the given feature dictionary.

    Args:
      feature (dict): The feature dictionary to add the baktfold key to.

    Returns:
      None

    Examples:
      >>> feature = {'sequence': 'ATG', 'start': 1, 'stop': 3, 'strand': '+'}
      >>> mark_as_baktfold(feature)
      >>> feature
      {'sequence': 'ATG', 'start': 1, 'stop': 3, 'strand': '+', 'baktfold': True}
    """
    # logger.info(
    #     f'baktfold found hit(s) for: seq={feature['sequence']}, start={feature['start']}, stop={feature['stop']}, strand={feature['strand']}'
    # )
    feature['baktfold'] = True


# def get_adjacent_genes(feature: dict, features: Sequence[dict], neighbors=3):
#     for idx, feat in enumerate(features):
#         if feat['locus'] == feature['locus']:
#             upstream_genes = []
#             if(idx >= 1):
#                 start = idx - neighbors
#                 if(start < 0 ):
#                     start = 0
#                 upstream_genes = features[start:idx]
#             downstream_genes = []
#             if(idx + 1 < len(features)):
#                 end = idx + 1 + neighbors
#                 if(end > len(features)):
#                     end = len(features)
#                 downstream_genes = features[idx+1:end]
#             upstream_genes.extend(downstream_genes)
#             for gene in upstream_genes:
#                 logger.debug(
#                     'extracted neighbor genes: seq=%s, start=%i, stop=%i, gene=%s, product=%s',
#                     gene['sequence'], gene['start'], gene['stop'], gene.get('gene', '-'), gene.get('product', '-')
#                 )
#             return upstream_genes
#     return []


# def select_gene_symbols(features: Sequence[dict]):
#     improved_genes = []
#     for feat in [f for f in features if len(f.get('genes', [])) > 1]:  # all CDS/sORF with multiple gene symbols
#         old_gene_symbol = feat['gene']
#         gene_symbol_prefixes = set([symbol[:3] for symbol in feat['genes'] if len(symbol) > 3])
#         if(len(gene_symbol_prefixes) == 1 ):  # multiple gene symbols of the same prefix: 
#             product_parts = feat.get('product', '').split()
#             for gene_symbol in feat['genes']:
#                 protein_symbol = gene_symbol[0].upper() + gene_symbol[1:]
#                 if(protein_symbol == product_parts[-1]):  # gene symbol is last part of product description which often is a specific gene/protein name
#                     if(gene_symbol != old_gene_symbol):
#                         feat['gene'] = gene_symbol
#                         logger.info(
#                             'gene product symbol selection: seq=%s, start=%i, stop=%i, new-gene=%s, old-gene=%s, genes=%s, product=%s',
#                             feat['sequence'], feat['start'], feat['stop'], gene_symbol, old_gene_symbol, ','.join(feat['genes']), feat.get('product', '-')
#                         )
#                         improved_genes.append(feat)
#         else:  # multiple gene symbols of varying prefixes are available, e.g. acrS, envR
#             logger.debug(
#                 'select gene symbol: seq=%s, start=%i, stop=%i, gene=%s, genes=%s, product=%s',
#                 feat['sequence'], feat['start'], feat['stop'], feat.get('gene', '-'), ','.join(feat['genes']), feat.get('product', '-')
#             )
#             adjacent_genes = get_adjacent_genes(feat, features, neighbors=3)
#             adjacent_gene_symbol_lists = [gene.get('genes', []) for gene in adjacent_genes]
#             adjacent_gene_symbols = [item for sublist in adjacent_gene_symbol_lists for item in sublist]  # flatten lists
#             adjacent_gene_symbol_prefixes = [gene_symbol[:3] for gene_symbol in adjacent_gene_symbols if len(gene_symbol) > 3]  # extract gene symbol prefixes, e.g. tra for traI, traX, traM
#             adjacent_gene_symbol_prefix_counts = {}
#             for gene_symbol_prefix in adjacent_gene_symbol_prefixes:
#                 if gene_symbol_prefix in adjacent_gene_symbol_prefix_counts:
#                     adjacent_gene_symbol_prefix_counts[gene_symbol_prefix] += 1
#                 else:
#                     adjacent_gene_symbol_prefix_counts[gene_symbol_prefix] = 1
#             logger.debug('neighbor gene symbol prefix counts: %s', adjacent_gene_symbol_prefix_counts)
#             count = 0
#             selected_gene_symbol = old_gene_symbol
#             for gene_symbol in feat['genes']:
#                 gene_symbol_prefix = gene_symbol[:3]
#                 gene_symbol_count = adjacent_gene_symbol_prefix_counts.get(gene_symbol_prefix, 0)
#                 if gene_symbol_count > count:  # select symbol if its prefix is dominant in the gene neighborhood (neihboorhood of 3 genes up-/downstream as operon proxy)
#                     selected_gene_symbol = gene_symbol
#                     count = gene_symbol_count
#             if(selected_gene_symbol != old_gene_symbol):
#                 feat['gene'] = selected_gene_symbol
#                 logger.info(
#                     'gene neighborhood symbol selection: seq=%s, start=%i, stop=%i, new-gene=%s, old-gene=%s, genes=%s, product=%s',
#                     feat['sequence'], feat['start'], feat['stop'], selected_gene_symbol, old_gene_symbol, ','.join(feat['genes']), feat.get('product', '-')
#                 )
#                 improved_genes.append(feat)
#     return improved_genes

###
# proteins
###


def annotate_aa(aas: Sequence[dict]):
    """
    Combines IPS and PSC annotations and marks hypotheticals.

    Args:
      aas (Sequence[dict]): A sequence of amino acid dictionaries to annotate.

    Returns:
      None

    Examples:
      >>> aas = [{'sequence': 'ATG', 'start': 1, 'stop': 3, 'strand': '+'}]
      >>> annotate_aa(aas)
      >>> aas
      [{'sequence': 'ATG', 'start': 1, 'stop': 3, 'strand': '+'}]
    """

    print('\tcombine annotations and mark hypotheticals...')

    for aa in aas:
        print(aa)
        combine_annotation(aa)  # combine IPS & PSC annotations and mark hypothetical
    logger.debug('analyze hypotheticals')
    hypotheticals = [aa for aa in aas if 'hypothetical' in aa]
    if(len(hypotheticals) > 0):
        print(f'\tanalyze hypothetical proteins: {len(hypotheticals)}')
        print('\tcalculated proteins statistics')
