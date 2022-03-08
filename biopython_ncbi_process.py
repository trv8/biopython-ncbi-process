import os
import shutil
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd

#Email (required) and API key (optional) used for NCBI Entrez queries
Entrez.email = 'email@example.com'
#Entrez.api_key = '0000'

#Path folders (2 required) for storing output XMLs: pmid_list XMLs and article_summaries XMLs
path_pmid = '../pmid'
path_articles = '../articles'

#Search term for generating PMIDs in export_xml_pmids_given_mesh
default_search_term = '"Alzheimer Disease"[Mesh]'

#Number of most recent articles to query
n_articles_max = 100

#Cache results
paramEutils = { 'usehistory':'N' }

################################
# Function Initialization
################################

def export_xml_pmids_given_searchterm(file_path = path_pmid, t=default_search_term, count=n_articles_max, batch_size=10000):
    #For a search string 'term', returns the latest  articles matching that term, outputs pmids to xml
    #Splits up requests in batches of size "increment"

    file_id = 0
    
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        r_max = batch_size
        num_left = count-start
        if num_left < batch_size:
            r_max = num_left
        
        handle = Entrez.esearch(retstart=start, retmax=r_max, term=t, sort='pub date', db='pubmed', retmode='xml',**paramEutils)
        record = handle.read()
        handle.close()
        file_name = file_path + '/pmid_list_' + str(file_id) + '.xml'
        
        if not os.path.isfile(file_name):
            file = open(file_name, 'wb')
            file.write(record)
            file.close()
            print("Successfully downloaded PMIDs %i to %i" % (start + 1, end))
        else:
            print("File already exists. Skipped file write for PMIDs %i to %i" % (start + 1, end))
        
        file_id += 1

        
def import_pmid_xml_batch(file_path = path_pmid):
    #outputs a list of PMIDs for XML files contained in a folder "file_path", which were generated from export_xml_pmids_given_searchterm
    files = os.listdir(file_path)
    if '.DS_Store' in files:
        files.remove('.DS_Store')
    
    PMIDs = []
    
    for f in files:
        f_name = file_path + '/' + f
        file = open(f_name, 'rb')
        record = Entrez.read(file)
        file.close()
        for p in record['IdList']:
            PMIDs.append(p)

    return PMIDs


def import_xml_entrez(file_name):
    #Returns Bio.Entrez object for processing for given input XML file at path 'path/file.xml'
    file = open(file_name, 'rb')
    record = Entrez.read(file)
    file.close()
    return record


def get_xml_articles_given_pmids(PMID_list=[], file_path = path_articles, count=n_articles_max, batch_size=500 ):
    #Given a python containing PMIDs, ['1234', '5678', ...], downloads articles as xml files to file_path_output, max # of 'count' articles

    def list_to_string(PMID_list):
        #Converts input PMIDs ['1234','5678'] to output string '1234,5678' for bioentrez parsing
        s = ''
        len_PMID = len(PMID_list)
        for i in range(len_PMID):
            s += PMID_list[i]
            if i< (len_PMID-1):
                s += ','
        return s

    file_id = 0
    
    for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            r_max = batch_size
            num_left = count-start
            if num_left < batch_size:
                r_max = num_left
                
            ids_to_fetch = ','.join(P[start:end])
            
            handle = Entrez.efetch(id=ids_to_fetch, db='pubmed', rettype='xml')
            record = handle.read()
            handle.close()
            file_name = file_path + '/pmid_articles_' + str(file_id) + '.xml'
            
            if not os.path.isfile(file_name):
                file = open(file_name, 'wb')
                file.write(record)
                file.close()
                print("Successfully downloaded articles %i to %i" % (start + 1, end))
            else:
                print("File already exists. Skipped file write for articles %i to %i" % (start + 1, end))
            
            file_id += 1

def merge_articles(file_path = path_articles):        
    files = os.listdir(file_path)
    
    if '.DS_Store' in files:
        files.remove('.DS_Store')
        
    if len(files) < 2:
        return 0

    #Copy first xml to merged file
    source_filename = file_path + '/' + files[0]
    merged_filename = file_path + '/' + 'merged.xml'
    shutil.copyfile(source_filename, merged_filename)

    #For additional xmls, append to merged
    for i in range(len(files)-1):
        tree_main = ET.parse(merged_filename)
        root_main = tree_main.getroot()
        
        tree = ET.parse(file_path + '/' + files[i+1])
        root = tree.getroot()

        for article in root:
            root_main.append(article)

        tree_main.write(merged_filename)
        print('Successfully merged ' + files[i+1])

    tree = ET.parse(merged_filename)
    root = tree.getroot()
    
    with open(merged_filename, "w", encoding='UTF-8') as xf:
        doc_type = '<?xml version="1.0" ?>\
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2019//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">'
        tostring = ET.tostring(root).decode('utf-8')
        file = f"{doc_type}{tostring}"
        xf.write(file)


################################
# End Function Initialization
################################

#Generate XML files for given search parameters
export_xml_pmids_given_searchterm()
P = import_pmid_xml_batch()
get_xml_articles_given_pmids(PMID_list=P)
merge_articles()

#Import records from merged.xml to test_record object
merged_filename = path_articles + '/' + 'merged.xml'
test_record = import_xml_entrez(merged_filename)
num_articles = len(test_record['PubmedArticle'])

#Prints the DescriptorNames for the articles, a:
for a in range(num_articles):
    PMID = str(test_record['PubmedArticle'][a]['MedlineCitation']['PMID'])
    print('\n'+'Article #' + str(a+1) + ' [PMID=' + PMID+']')
    mesh_article = test_record['PubmedArticle'][a]['MedlineCitation']['MeshHeadingList']
    for i in range(len(mesh_article)):
        print(str(mesh_article[i]['DescriptorName']))


#Access mesh ID instead of string descriptor
#dict(record2['PubmedArticle'][0]['MedlineCitation']['MeshHeadingList'][0]['QualifierName'][0].attributes.items())['UI']
