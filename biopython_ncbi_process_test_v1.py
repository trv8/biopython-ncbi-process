from Bio import Entrez
import os

#Email (required) and API key (optional) used for NCBI Entrez queries
Entrez.email = 'youremail@example.com'
#Entrez.api_key = 'put key here'

#Paths for storing output XMLs: PMID_list.xml and article_summaries.xml
path_pmid_output_given_mesh = 'choose path on HD /pmid_list.xml'
path_articles_given_pmid = 'choose path on HD /test_asummaries.xml'

#Search term for generating PMIDs in export_xml_pmids_given_mesh
default_search_term = '"Alzheimer Disease"[Mesh]'

#Number of most recent articles to query
n_articles_max = 5

#Cache results
paramEutils = { 'usehistory':'Y' }

################################
# Function Initialization
################################

def export_xml_pmids_given_mesh(file_name = path_pmid_output_given_mesh, t=default_search_term, n=n_articles_max):
    #For a search string 'term', returns the latest  articles matching that term, outputs pmids to xml
    handle = Entrez.esearch(retmax=n, term=t, sort='pub date', db='pubmed', retmode='xml',**paramEutils)
    record = handle.read()
    handle.close()

    if not os.path.isfile(file_name):
        file = open(file_name, 'wb')
        file.write(record)
        file.close()
    else:
        print('File already exists.')
    print('Done.')


def xml_pmids_to_list_object(file_name = path_pmid_output_given_mesh):
    #Imports an xml file from export_xml_pmids_given_mesh, outputs a python list containing the PMIDs
    temp = import_xml_entrez(file_name)
    return temp['IdList']


def get_xml_articles_given_pmids(PMID_list=[], file_name = path_articles_given_pmid):
    #Given a python containing PMIDs, ['1234', '5678', ...], downloads articles to xml file

    def list_to_string(PMID_list):
        #Converts input PMIDs ['1234','5678'] to output string '1234,5678' for bioentrez parsing
        s = ''
        len_PMID = len(PMID_list)
        for i in range(len_PMID):
            s += PMID_list[i]
            if i< (len_PMID-1):
                s += ','
        return s
    
    handle = Entrez.efetch(db='pubmed', id=list_to_string(PMID_list), sort='pub date', retmode='xml',**paramEutils)
    record = handle.read()
    handle.close()

    if not os.path.isfile(file_name):
        file = open(file_name, 'wb')
        file.write(record)
        file.close()
    else:
        print('File already exists.')
    print('Done.')


def import_xml_entrez(file_name):
    #Returns Bio.Entrez object for processing for given input XML file at file_name path /.xml
    file = open(file_name, 'rb')
    record = Entrez.read(file)
    file.close()
    return record


def list_mesh(file_name = path_articles_given_pmid):
    #Extracts the mesh components of article XMLs, to python list
    record = import_xml_entrez(file_name)
    mesh = record['PubmedArticle'][0]['MedlineCitation']['MeshHeadingList']
    return mesh


################################
# End Function Initialization
################################

#Ccomment the following 2 lines after generating xml files:
export_xml_pmids_given_mesh()
get_xml_articles_given_pmids(xml_pmids_to_list_object())

#Prints mesh for first article
test_record = import_xml_entrez(path_articles_given_pmid)
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








