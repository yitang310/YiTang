#!/usr/bin/env python
# coding: utf-8

import os
import json
import glob
import random
import numpy as np
import pandas as pd
from nltk import word_tokenize,pos_tag
from nltk.corpus import stopwords,wordnet
from nltk.stem import PorterStemmer,WordNetLemmatizer
from sentistrength import PySentiStr
os.chdir("App reviews procession and classification")

data=pd.read_csv("data3.csv")

len(data.index)

random.seed(553)
sampleindex=random.sample(range(len(data.index)),384)
data384=data.iloc[sampleindex]
data384.to_csv("data384.csv",index=None)

data1=pd.read_csv("output/Dataset3.csv",sep=",",encoding='unicode_escape')

#Remove the non-ASCII characters
data1["processed_text2"]=data1["processed_text"].str.encode("ascii", "ignore").str.decode("ascii")
#lower case
data1["processed_text2"]=data1["processed_text2"].apply(lambda x: x.lower())

data1.head(5)

#sentiScore total
def sentistr(x):
    senti = PySentiStr()
    senti.setSentiStrengthPath("SentiStrength.jar")
    senti.setSentiStrengthLanguageFolderPath("SentStrength_Data")
    result = senti.getSentiment(x,score='trinary') #positive rating, negative rating and neutral rating
    return result
sentis=data1["processed_text2"].copy()
sentis=sentis.apply(sentistr)
def senti_score(x):
    if x[0][2]>=0:
        return x[0][0]
    else:
        return x[0][1]


data_total=pd.DataFrame()

#comment
data_total["comment"]=data1["text"]

#rating
data_total["rating"]=data1["score"]

#past
def tense_past(x):
    text = word_tokenize(x)
    tagged = pos_tag(text)
    return len([word for word in tagged if word[1] in ["VBD", "VBN"]])
data_total["past"]=data1["processed_text2"].apply(tense_past)

#stopwords_removal
CUSTOM_STOPWORDS = ['i', 'me','up','my', 'myself', 'we', 'our', 'ours',
                    'ourselves', 'you', 'your', 'yours','yourself', 'yourselves',
                    'he', 'him', 'his', 'himself', 'she', 'her', 'hers' ,'herself',
                    'it', 'its', 'itself', 'they', 'them', 'their', 'theirs',
                    'themselves' ,'am', 'is', 'are','a', 'an', 'the', 'and','in',
                    'out', 'on','up','down', 's', 't']
def remove_stopwords(x):
    sw_cs=" ".join([word for word in x.split() if word not in (CUSTOM_STOPWORDS)])
    return sw_cs
data_total["stopwords_removal"]=data1["processed_text2"].apply(remove_stopwords)

#reviewer
data_total["reviewer"]=None

#id
data_total["id"]=data1["id"]

#stemmed
def stemming(x):
    ps = PorterStemmer()
    s_s=" ".join([ps.stem(word) for word in x.split()])
    return s_s
data_total["stemmed"]=data1["processed_text2"].apply(stemming)

#fee
data_total["fee"]=None

#title
data_total["title"]=None

#label
##1. bug reports 2.feature requests 3. user experiences 4.text ratings
def class1234(x):
    x=int(x)
    if x==1:
        return "Bug"
    elif x==2:
        return "Feature"
    elif x==3:
        return "UserExperience"
    elif x==4:
        return "Rating"
data_total["label"]=data1["final"].apply(class1234)

#future
def tense_future(x):
    text = word_tokenize(x)
    tagged = pos_tag(text)
    return len([word for word in tagged if word[1] == "MD"])
data_total["future"]=data1["processed_text2"].apply(tense_future)

#lemmatized_comment
def nltk2wn_tag(nltk_tag):
    if nltk_tag.startswith('J'):
        return wordnet.ADJ
    elif nltk_tag.startswith('V'):
        return wordnet.VERB
    elif nltk_tag.startswith('N'):
        return wordnet.NOUN
    elif nltk_tag.startswith('R'):
        return wordnet.ADV
    else:
        return None
def lemmatize_sentence(x):
    lemmatizer = WordNetLemmatizer()
    nltk_tagged = pos_tag(word_tokenize(x))
    text = word_tokenize(x)
    tagged = pos_tag(text)
    wn_tagged = map(lambda y: (y[0], nltk2wn_tag(y[1])), tagged)
    res_words = []
    for word, tag in wn_tagged:
        if tag is None:
            res_words.append(word)
        else:
            res_words.append(lemmatizer.lemmatize(word, tag))
    return " ".join(res_words)
data_total["lemmatized_comment"]=data1["processed_text2"].apply(lemmatize_sentence)

#sentiScore
data_total["sentiScore"]=sentis.apply(senti_score)

#sentiScore_neg
data_total["sentiScore_neg"]=sentis.apply(lambda x : x[0][1])

#reviewId
data_total["reviewId"]=data1["id"]

#stopwords_removal_nltk
def remove_stopwords_nltk(x):
    sw_s=" ".join([word for word in x.split() if word not in (stopwords.words('english'))])
    return sw_s
data_total["stopwords_removal_nltk"]=data1["processed_text2"].apply(remove_stopwords_nltk)

#present_simple
def tense_presentsimple(x):
    text = word_tokenize(x)
    tagged = pos_tag(text)
    return len([word for word in tagged if word[1] in ["VBP", "VBZ"]])
data_total["present_simple"]=data1["processed_text2"].apply(tense_presentsimple)

#dataSource
data_total["dataSource"]=data1["datetime"] + "_TOP_FREE"+data1["category"]+"_APPS"

#appId
data_total["appId"]=data1["appTitle"]

#date
data_total["date"]=data1["date"]

#sentiScore_pos
data_total["sentiScore_pos"]=sentis.apply(lambda x : x[0][0])

#present_con
def tense_presentcon(x):
    text = word_tokenize(x)
    tagged = pos_tag(text)
    return len([word for word in tagged if word[1] == "VBG"])
data_total["present_con"]=data1["processed_text2"].apply(tense_presentcon)

#length_words
data_total["length_words"]=data1["len"]

#stopwords_removal_lemmatization
data_total["stopwords_removal_lemmatization"]=data_total["stopwords_removal"].apply(lemmatize_sentence)


# In[55]:


data_total.groupby(["label"]).size()


# In[61]:


#seperate data
random.seed(553)
arr4=np.array(["Bug","Feature","UserExperience","Rating"])
for i in arr4:
    tmp4=data_total[data_total["label"]==i]
    tmp4not=data_total[data_total["label"]!=i]
    if len(tmp4)<=len(tmp4not):
        index=random.sample(range(len(tmp4not.index)),len(tmp4.index))
        tmp4not=tmp4not.iloc[index]
        tmp4not["label"]=tmp4not["label"].apply(lambda x: "Not_"+i)
        tmp4_total=pd.concat([tmp4,tmp4not],ignore_index=True)
    else:
        index=random.sample(range(len(tmp4.index)),len(tmp4not.index))
        tmp4=tmp4.iloc[index]
        tmp4not["label"]=tmp4not["label"].apply(lambda x: "Not_"+i)
        tmp4_total=pd.concat([tmp4,tmp4not],ignore_index=True)
    print(i,len(tmp4_total.index))
    #to json
    data_dict=tmp4_total.to_dict("records")
    with open("output/data/"+i+"_tt_own.json","w") as f:
        json.dump(data_dict,f)



data_total.head()


#to json
data_dict=data_total.to_dict("records")
with open("output/data/all_own.json", "w") as f:
    json.dump(data_dict,f)



#merge json
arr4=np.array(["Bug","Feature","UserExperience","Rating"])
for i in arr4:
    result=[]
    for f in glob.glob("**/*"+i+"*.json", recursive=True):
        result.append(f)
    with open(result[0], "r") as infile1:
        re_l1=json.load(infile1)
    with open(result[1], "r") as infile2:
        re_l2=json.load(infile2)
    re_l=re_l1+re_l2
    with open("data/"+i+"_merged.json", "w") as outfile:
        json.dump(re_l, outfile)
