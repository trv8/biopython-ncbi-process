import os
import shutil
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd

#Email (required) and API key (optional) used for NCBI Entrez queries
Entrez.email = 'email@example.com'
Entrez.api_key = '0000'

#Path folders (2 required) for storing output XMLs: pmid_list XMLs and article_summaries XMLs
path_pmid = './entrez_data/pmids'
path_articles = './entrez_data/articles'

#Search term for generating PMIDs in export_xml_pmids_given_mesh
default_search_term = 'Alzheimer'

#Number of most recent articles to query
n_articles_max = 25000

#Exclude articles that are not MESH tagged?
exclude_untagged = True

#Cache results
paramEutils = { 'usehistory':'Y' }


#XML namespace (required for XML merge function)
ET.register_namespace('mml', "http://www.w3.org/1998/Math/MathML")

################################
# Function Initialization
################################

def export_xml_pmids_given_searchterm(file_path = path_pmid, t=default_search_term, count=n_articles_max, batch_size=10000):
    #For a search string 'term', returns the latest PMIDs matching that term, outputs pmids to xml
    #Splits up requests in batches of size "batch_size"

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


def get_xml_articles_given_pmids(PMID_list=[], file_path = path_articles, count=n_articles_max, batch_size=1000, **paramEutils ):
    #Given a python list containing PMIDs, ['1234', '5678', ...], downloads articles as xml files to file_path_output, max # of 'count' articles

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

def merge_articles(file_path = path_articles, output_file = 'merged.xml'):
    #Merges XML files contained in file_path folder to an output file named 'merged.xml'
    files = os.listdir(file_path)
    
    if '.DS_Store' in files:
        files.remove('.DS_Store')
        
    if len(files) < 2:
        return 0

    #Copy first xml to merged file
    source_filename = file_path + '/' + files[0]
    merged_filename = file_path + '/' + output_file
    shutil.copyfile(source_filename, merged_filename)

    #For additional xmls, append to merged
    for i in range(len(files)-1):
        tree_main = ET.parse(merged_filename)
        root_main = tree_main.getroot()
        
        tree = ET.parse(file_path + '/' + files[i+1])
        root = tree.getroot()

        for article in root:
            root_main.append(article)

        tree_main.write(merged_filename, method = 'xml',encoding = 'utf-8')
        print('Successfully merged ' + files[i+1])



    tree = ET.parse(merged_filename)
    root = tree.getroot()
    
    with open(merged_filename, "w", encoding='UTF-8') as xf:
        doc_type = '<?xml version="1.0" ?>\
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2019//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">'
        tostring = ET.tostring(root).decode('utf-8')
        file = f"{doc_type}{tostring}"
        xf.write(file)
        


def pd_dataframe_xml_multiple(file_path = path_articles, exclude_utag = exclude_untagged):
    #Returns a Pandas dataframe containing relevant information from XML files located at file_path folder
    
    files = os.listdir(file_path)
    if '.DS_Store' in files:
        files.remove('.DS_Store')
    if 'merged.xml' in files:
        files.remove('merged.xml')

    PMID_df = []
    Mesh_df = []
    Desc_names_df = []
    All_tags_df = []
    All_tags_no_major_df = []
        
    for file in files:
        file_dir = file_path + '/' + file
        record = import_xml_entrez(file_dir)
        num_articles = len(record['PubmedArticle'])

        for a in range(num_articles):
            PMID = str(record['PubmedArticle'][a]['MedlineCitation']['PMID'])
            is_tagged = 'MeshHeadingList' in record['PubmedArticle'][a]['MedlineCitation'].keys()

            if is_tagged:
                mesh_of_article = record['PubmedArticle'][a]['MedlineCitation']['MeshHeadingList']
            else:
                mesh_of_article = []
            
            mesh_temp = []
            desc_names_temp = []
            all_tags_temp = []
            all_tags_no_major_temp = []
            
            for i in range(len(mesh_of_article)):
                mesh_temp.append(mesh_of_article[i])
                desc_names_temp.append((str(mesh_of_article[i]['DescriptorName'])))


                #If there are qualifiers present
                if len(mesh_of_article[i]['QualifierName']) != 0:
                    for j in range(len(mesh_of_article[i]['QualifierName'])):
                        tag_str = str(mesh_of_article[i]['DescriptorName']) + '/' + mesh_of_article[i]['QualifierName'][j]
                        tag_str_no_major = tag_str
                        
                        #Add * to subtopic if it's a major topic
                        if dict(mesh_of_article[i]['QualifierName'][j].attributes.items())['MajorTopicYN'] == 'Y':
                            tag_str += '*' 
                        
                        all_tags_temp.append(tag_str)
                        all_tags_no_major_temp.append(tag_str_no_major)
                else:
                    tag_str = str(mesh_of_article[i]['DescriptorName'])
                    tag_str_no_major = tag_str

                    #Add * to main descriptor if it's a major topic itself
                    if dict(mesh_of_article[i]['DescriptorName'].attributes.items())['MajorTopicYN'] == 'Y':
                        tag_str += '*'
                        
                    all_tags_temp.append(tag_str)
                    all_tags_no_major_temp.append(tag_str_no_major)
                    
            if not(not(is_tagged) and exclude_utag):
                PMID_df.append(PMID)                   
                Mesh_df.append(mesh_temp)
                Desc_names_df.append(desc_names_temp)
                All_tags_df.append(all_tags_temp)
                All_tags_no_major_df.append(all_tags_no_major_temp)
                
    d = {'pmid': PMID_df, 'mesh': Mesh_df, 'major_tags': Desc_names_df, 'all_tags': All_tags_df,'all_tags_no_major': All_tags_no_major_df}
    df = pd.DataFrame(data=d)
    
    return df

def pd_freq_table(df, column = 'major_tags'):
    #Returns a Pandas data frame containing Frequency data for tags in tag columns all_tags, major_tags, etc.
    num_articles = len(df)
    Temp_list = []
    for i in range(len(df[column])):
        for j in range(len(df[column][i])):
            Temp_list.append(df[column][i][j])

    Temp_list = pd.Series(Temp_list)
    Freq_df = Temp_list.value_counts()
    Freq_df = pd.DataFrame({column:Freq_df.index, 'count':Freq_df.values})
    Freq_df['freq'] = Freq_df['count']/num_articles
    
    return Freq_df

################################
# End Function Initialization
################################

#Generate XML files for given search parameters
##export_xml_pmids_given_searchterm()
##P = import_pmid_xml_batch()
##get_xml_articles_given_pmids(PMID_list=P)

#Use to generate main dataframe. Can comment and import .pkl files for faster use after generation
#Using df.to_json('df.json') and pd.read_json('df.json')

#df = pd_dataframe_xml_multiple()
#df.to_json('df.json')

df = pd.read_json('df.json')
df_freq = pd_freq_table(df, column = 'major_tags')

print(df_freq)
