import nltk
from bs4 import BeautifulSoup
import re
import nltk
from nltk.corpus import stopwords
from sklearn.feature_extraction.text import CountVectorizer

filename="C:/Users/halidziya/Desktop/word2vec/realtimewwii.txt"





with open(filename,'r',encoding='utf-8') as infile:
    text = infile.read();
    text = re.sub("[^a-zA-Z0-9\n]"," ",text).lower() #Clean
    text = text.split("\n")
    cleaned = ""
    i=0

    for aline in text: #Stop words
        words = aline.split()  
        #words = [w for w in words if not w in stopwords.words("english")] 
        clean = " ".join( words )
        text[i] = clean
        i=i+1
    vectorizer = CountVectorizer(analyzer = "word",    tokenizer = None,  min_df=5,  preprocessor = None, stop_words = None, max_features = 5000) 
    train_data_features = vectorizer.fit_transform(text)
    train_data_features = train_data_features.toarray()
    vocab = vectorizer.get_feature_names()


import csv
with open("bow.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(train_data_features)

import csv
with open("vocabulary.txt", "w") as f:
    writer = csv.writer(f)
    writer.writerows(vocab)

with open("cleaned.txt","w") as f:
   for item in text:
       f.writelines('%s\n' % item)