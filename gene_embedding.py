#import Python libaries needed for training embbeded vectors 
from gensim.models.word2vec import LineSentence
from gensim.models import Word2Vec
import time # for checking how long the training process takes

# where the data is located
input_data = 'pubmed/pubmed.txt'

# parameters for training
sg_ = 1 # the training algorithm. If sg=0, CBOW is used. Otherwise (sg=1), skip-gram is employed.
size_ = 50 #  the dimensionality of the feature vectors
window_ = 5 # the context size or the maximum distance between the current and predicted word

# keep the time starting the training
start_time = time.time()

# training embedded vectors for the dataset with the parameters specified above
model = Word2Vec(LineSentence(input_data), sg = sg_, size = size_, window = window_)
# save the model learned into model file
model.save('pubmed/gene_embedding.embed')

# show how long does it take to train the word vectors
runtime = time.time() - start_time
print("--- Running time: %s seconds ---" % (runtime))