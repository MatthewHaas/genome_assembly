punctuation_chars = ["'", '"', ",", ".", "!", ":", ";", '#', '@']
# lists of words to use
positive_words = []
with open("positive_words.txt") as pos_f:
    for lin in pos_f:
        if lin[0] != ';' and lin[0] != '\n':
            positive_words.append(lin.strip())


negative_words = []
with open("negative_words.txt") as pos_f:
    for lin in pos_f:
        if lin[0] != ';' and lin[0] != '\n':
            negative_words.append(lin.strip())

def strip_punctuation(string):
    for char in string:
        if char in punctuation_chars:
            string = string.replace(char, "")
    return string

def get_pos(string):
    pos_words = 0
    word_list = string.split()
    for word in word_list:
        word = word.lower()
        word = strip_punctuation(word)
        if word in positive_words:
            pos_words += 1
        else:
            continue
    return pos_words

def get_neg(string):
    neg_words = 0
    word_list = string.split()
    for word in word_list:
        word = word.lower()
        word = strip_punctuation(word)
        if word in negative_words:
            neg_words += 1
        else:
            continue
    return neg_words

raw_data = open('project_twitter_data.csv','r')

with open('resulting_data.csv', 'w') as outfile:
    outfile.write("Number of Retweets, Number of Replies, Positive Score, Negative Score, Net Score")
    outfile.write("\n")
    infile = raw_data.readlines()
    header = infile.pop(0)
    for row in infile:
        line = row.strip().split(',')
        outfile.write('{}, {}, {}, {}, {}'.format(line[1], line[2], get_pos(line[0]), get_neg(line[0]), (get_pos(line[0])-get_neg(line[0]))))  
        outfile.write('\n')

# Close files        
raw_data.close()
outfile.close()