import os
import shutil
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
from itertools import combinations
import time

#Email (required) and API key (optional) used for NCBI Entrez queries
Entrez.email = 'email@example.com'
Entrez.api_key = '000000'

#Path folders (2 required) for storing output XMLs: pmid_list XMLs and article_summaries XMLs
path_pmid = './entrez_data/pmids'
path_articles = './entrez_data/articles'

#Search term for generating PMIDs in export_xml_pmids_given_mesh
default_search_term = 'Alzheimer'

#Number of most recent articles to query
n_articles_max = 194750

#Exclude articles that are not MESH tagged?
#SET TO FALSE FOR NOW. DO NOT CHANGE VALUE
exclude_untagged = False

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
        
        handle = Entrez.esearch(db='pubmed',retstart=start, retmax=r_max, term=t, sort='pub date', retmode='xml',**paramEutils)
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
                
            ids_to_fetch = ','.join(PMID_list[start:end])
            
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
        


def pd_dataframe_xml_multiple(file_path = path_articles, exclude_utag = exclude_untagged, gen_major_qualifiers = False, gen_append_mesh = False):
    #Returns a Pandas dataframe containing relevant information from XML files located at file_path folder
    #Setting gen_major_qualifiers = True will generate an additional column which appends major topics with * (i.e. Smoking*)
    #Setting gen_append_mesh = True will generate an additional column which appends the unfiltered mesh data for custom processing
    
    files = os.listdir(file_path)
    if '.DS_Store' in files:
        files.remove('.DS_Store')
    if 'merged.xml' in files:
        files.remove('merged.xml')

    PMID_df = []
    Mesh_df = []
    Desc_names_df = []
    All_tags_no_major_df = []
    All_tags_df = []
    Keywords_df = []
    RN_df = []
        
    for file in files:
        print('Importing data from ' + file + '...')
        
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
            kw_temp = []
            rn_temp = []

            #If there are keywords present
            if len(record['PubmedArticle'][a]['MedlineCitation']['KeywordList']) != 0:
                #for every keyword
                for i in range(len(record['PubmedArticle'][a]['MedlineCitation']['KeywordList'][0])):
                    #add keyword string to temporary list
                    kw_temp.append(str(record['PubmedArticle'][a]['MedlineCitation']['KeywordList'][0][i]))

            #if there are RN values present
            if 'ChemicalList' in record['PubmedArticle'][a]['MedlineCitation'].keys():
                
                RN_of_a = record['PubmedArticle'][a]['MedlineCitation']['ChemicalList']
                
                #If the RN field isn't empty
                if len(RN_of_a) != 0:
                    for i in range(len(RN_of_a)):
                        rn_name_string = str(RN_of_a[i]['NameOfSubstance'])
                        rn_id_string = str(RN_of_a[i]['RegistryNumber'])
                        rn_str_temp = 'RN: ' + rn_id_string + ' (' + rn_name_string + ')'
                        rn_temp.append(rn_name_string)

                        
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
                        
                        if gen_major_qualifiers:
                            all_tags_temp.append(tag_str)
                            
                        all_tags_no_major_temp.append(tag_str_no_major)
                else:
                    tag_str = str(mesh_of_article[i]['DescriptorName'])
                    tag_str_no_major = tag_str

                    #Add * to main descriptor if it's a major topic itself
                    if dict(mesh_of_article[i]['DescriptorName'].attributes.items())['MajorTopicYN'] == 'Y':
                        tag_str += '*'
                        
                    if gen_major_qualifiers:
                        all_tags_temp.append(tag_str)
                        
                    all_tags_no_major_temp.append(tag_str_no_major)
                    
            if not(exclude_utag):
                PMID_df.append(PMID)                   
                Desc_names_df.append(desc_names_temp)
                All_tags_no_major_df.append(all_tags_no_major_temp)
                Keywords_df.append(kw_temp)
                RN_df.append(rn_temp)
                
                if gen_append_mesh:
                    Mesh_df.append(mesh_temp)
                if gen_major_qualifiers:
                    All_tags_df.append(all_tags_temp)
                    
    d = {'pmid': PMID_df, 'main_tags': Desc_names_df,'all_tags_no_major_topic': All_tags_no_major_df, 'keywords': Keywords_df, 'rn_values': RN_df}
    
    if gen_major_qualifiers:
        d['all_tags_with_qualifiers'] = All_tags_df
    if gen_append_mesh:
        d['mesh'] = Mesh_df
        
    df = pd.DataFrame(data=d)
    return df

def pd_freq_table(df, column = 'main_tags'):
    #Returns a Pandas data frame containing Frequency data for lists of tags in tag columns all_tags, major_tags, etc.
    num_articles = len(df)
    Temp_list = []
    for i in range(len(df[column])):
        for j in range(len(df[column][i])):
            Temp_list.append(df[column][i][j])

    Temp_list = pd.Series(Temp_list)
    df_freq = Temp_list.value_counts()
    df_freq = pd.DataFrame({column:df_freq.index, 'count':df_freq.values})
    df_freq['freq'] = df_freq['count']/num_articles
    
    return df_freq

def pd_freq_table_remove_values(df_freq, f_upper_limit = 0.333333333333333, count_lower_limit = 5):
    #removes values from dataframe with freq > f_limit and count <= count_lower_limit
    df_freq_updated = df_freq
    df_freq_updated = df_freq_updated[df_freq_updated['freq'] < f_upper_limit]
    df_freq_updated = df_freq_updated[df_freq_updated['count'] > count_lower_limit]
    
    return df_freq_updated


def generate_pair_pmid(df):
    #returns pandas dataframe containing all pairs of PMIDs under column 'pmid'
    #must also contain column with explicit indices 'index'
    PMID_list = df['pmid']
    index_list = df['index']
    left_pmid = []
    right_pmid = []
    left_i = []
    right_i = []
    comb = combinations(range(len(PMID_list)),2)
    
    for c in comb:
        left_index = c[0]
        right_index = c[1]
        left_pmid.append(PMID_list[left_index])
        right_pmid.append(PMID_list[right_index])
        left_i.append(index_list[left_index])
        right_i.append(index_list[right_index])

    d = {'leftpmid': left_pmid, 'rightpmid': right_pmid, 'lefti': left_i, 'righti': right_i}
    df_pairs = pd.DataFrame(data=d)
    return df_pairs


def compute_distance_pmid_pair(df, df_pairs, df_freq):
    #Given a dataframe df_pairs containing pairs of pmids (leftpmid, rightpmid, lefti, righti)
    #and a main dataframe with R/L indices, computes distance metric
    scores = []

    #Generate an explicit dictionary of term <-> frequency pairs. Keys = term, values = freq
    F_dict = dict(zip(df_freq.merge_tag_rn, df_freq.freq))
    
    for p in range(len(df_pairs)):
        lefti_df = df_pairs['lefti'][p]
        righti_df = df_pairs['righti'][p]
        left_set = df.iloc[lefti_df]['merge_tag_rn']
        right_set = df.iloc[righti_df]['merge_tag_rn']

        shared_tags_list = list(set(left_set) & set(right_set))
        similarity_score = 0
        for i in range(len(shared_tags_list)):
            term = shared_tags_list[i]

            #This if statement allows for processing post-filtered sets
            if term in F_dict.keys():
                similarity_score += 1.0/F_dict[term]

        similarity_score = round(similarity_score)
        scores.append(similarity_score)


    df_pairs_scored = df_pairs.assign(sim_score = scores)
    
    return df_pairs_scored
        

################################
# End Function Initialization
################################

#Generate XML files for given search parameters
#export_xml_pmids_given_searchterm()
#P = import_pmid_xml_batch()
#get_xml_articles_given_pmids(PMID_list=P)

#Use to generate main dataframe. Can comment and import .json files for faster use after generation
#Using df.to_json('df.json') and pd.read_json('df.json')

#df = pd_dataframe_xml_multiple()

#df.to_json('df.json')

df = pd.read_json('df.json')


d = []
df.insert(0,'merge_tag_rn','')
for i in range(len(df)):
    t = df['main_tags'][i]
    r = df['rn_values'][i]
    v = list(set(t+r))
    d.append(v)
    
df['merge_tag_rn'] = d



df2 = df.sample(n=100)
df2 = df2.reset_index()
print('Generating PMID pairs...')
df_pairs = generate_pair_pmid(df2)

print('Generating frequency table...')
df_freq = pd_freq_table(df, column = 'merge_tag_rn')

#Filters data
print('Filtering frequency data...')
df_freq_updated = pd_freq_table_remove_values(df_freq)

print('Scoring...')
df_pairs_scored = compute_distance_pmid_pair(df,df_pairs,df_freq_updated)



##import matplotlib.pyplot as plt
##df_pairs_scored['sim_score'].plot.hist(bins=1000)
##plt.show()


#df_freq = pd_freq_table(df, column = 'main_tags')

#print(df_freq)

#Remove items with too high frequency, f_limit, are removed and placed in df_freq_updated
#df_freq_updated = pd_freq_table_remove_values(df_freq)

#print(df_freq_updated)
