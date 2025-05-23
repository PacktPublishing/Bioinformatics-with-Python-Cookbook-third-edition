<p align="center"><a href="https://packt.link/mlsumgh"><img src="https://static.packt-cdn.com/assets/images/ML Summit Banner v3 1200x627.png" alt="Machine Learning Summit 2025"/></a></p>

## Machine Learning Summit 2025
**Bridging Theory and Practice: ML Solutions for Today’s Challenges**

3 days, 20+ experts, and 25+ tech sessions and talks covering critical aspects of:
- **Agentic and Generative AI**
- **Applied Machine Learning in the Real World**
- **ML Engineering and Optimization**

👉 [Book your ticket now >>](https://packt.link/mlsumgh)

---

## Join Our Newsletters 📬

### DataPro  
*The future of AI is unfolding. Don’t fall behind.*

<p><a href="https://landing.packtpub.com/subscribe-datapronewsletter/?link_from_packtlink=yes"><img src="https://static.packt-cdn.com/assets/images/DataPro NL QR Code.png" alt="DataPro QR" width="150"/></a></p>

Stay ahead with [**DataPro**](https://landing.packtpub.com/subscribe-datapronewsletter/?link_from_packtlink=yes), the free weekly newsletter for data scientists, AI/ML researchers, and data engineers.  
From trending tools like **PyTorch**, **scikit-learn**, **XGBoost**, and **BentoML** to hands-on insights on **database optimization** and real-world **ML workflows**, you’ll get what matters, fast.

> Stay sharp with [DataPro](https://landing.packtpub.com/subscribe-datapronewsletter/?link_from_packtlink=yes). Join **115K+ data professionals** who never miss a beat.

---

### BIPro  
*Business runs on data. Make sure yours tells the right story.*

<p><a href="https://landing.packtpub.com/subscribe-bipro-newsletter/?link_from_packtlink=yes"><img src="https://static.packt-cdn.com/assets/images/BIPro NL QR Code.png" alt="BIPro QR" width="150"/></a></p>

[**BIPro**](https://landing.packtpub.com/subscribe-bipro-newsletter/?link_from_packtlink=yes) is your free weekly newsletter for BI professionals, analysts, and data leaders.  
Get practical tips on **dashboarding**, **data visualization**, and **analytics strategy** with tools like **Power BI**, **Tableau**, **Looker**, **SQL**, and **dbt**.

> Get smarter with [BIPro](https://landing.packtpub.com/subscribe-bipro-newsletter/?link_from_packtlink=yes). Trusted by **35K+ BI professionals**, see what you’re missing.

# Important note

As of September 2023, sgkit used in Chapter 6 has a bug. If you use Chapter 6, `pip install` this:

`pip install "sgkit[plink] @ git+https://github.com/pystatgen/sgkit@72ba77697fe1ca975dc5761a0a9c4355e5e49419"`



# Bioinformatics-with-Python-Cookbook-third-edition

<a href="https://www.packtpub.com/product/bioinformatics-with-python-cookbook-third-edition/9781803236421"><img src="https://static.packt-cdn.com/products/9781803236421/cover/smaller" alt="Bioinformatics with Python Cookbook - Third Edition" height="256px" align="right"></a>

This is the code repository for [Bioinformatics with Python Cookbook - Third Edition](https://www.packtpub.com/product/bioinformatics-with-python-cookbook-third-edition/9781803236421), published by Packt.

**Use modern Python libraries and applications to solve real-world computational biology problems**

## What is this book about?
Bioinformatics is an active research field that uses a range of simple-to-advanced computations to extract valuable information from biological data, and this book will show you how to manage these tasks using Python.

This updated third edition of the Bioinformatics with Python Cookbook begins with a quick overview of the various tools and libraries in the Python ecosystem that will help you convert, analyze, and visualize biological datasets. Next, you'll cover key techniques for next-generation sequencing, single-cell analysis, genomics, metagenomics, population genetics, phylogenetics, and proteomics with the help of real-world examples. You'll learn how to work with important pipeline systems, such as Galaxy servers and Snakemake, and understand the various modules in Python for functional and asynchronous programming. This book will also help you explore topics such as SNP discovery using statistical approaches under high-performance computing frameworks, including Dask and Spark. In addition to this, you’ll explore the application of machine learning algorithms in bioinformatics.

By the end of this bioinformatics Python book, you'll be equipped with the knowledge you need to implement the latest programming techniques and frameworks, empowering you to deal with bioinformatics data on every scale.

This book covers the following exciting features: 
* Become well-versed with data processing libraries such as NumPy, pandas, arrow, and zarr in the context of bioinformatic analysis
* Interact with genomic databases
* Solve real-world problems in the fields of population genetics, phylogenetics, and proteomics
* Build bioinformatics pipelines using a Galaxy server and Snakemake
* Work with functools and itertools for functional programming
* Perform parallel processing with Dask on biological data
* Explore principal component analysis (PCA) techniques with scikit-learn

If you feel this book is for you, get your [copy](https://www.amazon.in/Bioinformatics-Python-Cookbook-bioinformatics-computational/dp/1789344697/ref=sr_1_2?keywords=Bioinformatics+with+Python+Cookbook+-+Third+Edition&qid=1665382032&sr=8-2) today!

<a href="https://www.packtpub.com/product/bioinformatics-with-python-cookbook-third-edition/9781803236421"><img src="https://raw.githubusercontent.com/PacktPublishing/GitHub/master/GitHub.png" alt="https://www.packtpub.com/" border="5" /></a>

## Instructions and Navigations
All of the code is organized into folders.

The code will look like the following:
```
from Bio import SeqIO
genome_name = 'PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta'
recs = SeqIO.parse(genome_name, 'fasta')
for rec in recs:
print(rec.description)
```
**Following is what you need for this book:**
This book is for bioinformatics analysts, data scientists, computational biologists, researchers, and Python developers who want to address intermediate-to-advanced biological and bioinformatics problems. Working knowledge of the Python programming language is expected. Basic knowledge of biology will also be helpful.

With the following software and hardware list you can run all code files present in the book (Chapter 1-12).

### Software and Hardware List

| Chapter  | Software required                                                                    | OS required                        |
| -------- | -------------------------------------------------------------------------------------| -----------------------------------|
|  	1-12	   | Python 3.9                             			  | Any OS | 		
|  	1-12	   | Numpy, Pandas and Matplotlib                             			  | Any OS | 		
|  	1-12	   | BioPython                             			  | Any OS | 		
|  	1-12	   | DAsk, Zarr, Sckit-learn                             			  | Any OS | 		

We also provide a PDF file that has color images of the screenshots/diagrams used in this book. [Click here to download it](https://packt.link/3KQQO).

  
## Get to Know the Author
**Tiago Antao** is a bioinformatician who is currently working in the field of genomics. A former computer scientist, Tiago moved into computational biology with an MSc in bioinformatics from the Faculty of Sciences at the University of Porto, Portugal, and a PhD on the spread of drug-resistant malaria from the Liverpool School of Tropical Medicine, UK. Post his doctoral, Tiago worked with human datasets at the University of Cambridge, UK and with mosquito whole-genome sequencing data at the University of Oxford, UK, before helping to set up the bioinformatics infrastructure at the University of Montana, USA. He currently works as a data engineer in the biotechnology field in Boston, MA. He is one of the co-authors of Biopython, a major bioinformatics package written in Python.
### Download a free PDF

 <i>If you have already purchased a print or Kindle version of this book, you can get a DRM-free PDF version at no cost.<br>Simply click on the link to claim your free PDF.</i>
<p align="center"> <a href="https://packt.link/free-ebook/9781803236421">https://packt.link/free-ebook/9781803236421 </a> </p>
