# App reviews procession and classification

### Dec. 2019 -- Feb. 2020

Goals: handling huge reviews which are mobile applications text with a total of 2715303 by natural language processing (NLP), and reproducing the program from report On the Automatic Classification of App Reviews with randomly selected App reviews to classify them into four categories.

Methods: Processing the reviews including removing non-English reviews, non-ASCII characters, punctuations, adjacent duplicate characters and reviews that contain two or fewer words by langid, string and regular expression for the original dataset, counting future tense words, past tense words, present simple tense words and present continuous tense words for each review by nltk and calculating sentiment rating for App reviews by sentistrength, selecting 384 rows App review text from the processed dataset with manually labelled four categories to test the accuracy

Tools: Manipulating Python to complete the project
