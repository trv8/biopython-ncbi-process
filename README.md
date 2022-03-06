# biopython-ncbi-process
Work in progress.

Currently takes most recent 5 articles from NCBI's PubMed database matching search term "Alzheimer Disease"[MeSH] and returns the associated MeSH tags with these articles. Stores XML data on local disk to minimize number of queries to the PubMed database.

Requirements:
Python 3.8 https://www.python.org/ ,
BioPython
https://biopython.org/docs/1.76/api/index.html#

How to use:

*Edit the following line to your email:
Entrez.email = 'your_email@example.com'

*Edit default_search_term to a search of choice.

*Edit path_pmid_output_given_mesh and path_articles_given_pmid to a location on your hard drive to store xml files.



Output:
MeSH tags associated with the most recent 5 articles are outputted to the terminal.
