## OBJECTIVES

Suppose that in a random sample of documents from the Web, the toponym "Tresjuncos" occurs together with Madrid, Valencia, and Barcelona 200 times, and with Mexico City, Merida, and Guadalajara only 20 times. Is Tresjuncos more likely to be in Spain or Mexico? Although this example is easy, the problem of locating toponyms based on their distributions in text is much more challenging in practice. 

GRAIL probabilistically estimates geographic locations from place name co-occurrences in text. It produces heatmap-like images which represent GRAIL's estimation of the location to which a particular toponym refers, based on the frequency with which that toponym co-occurs with others in the specified corpus. Areas colored in dark blue have a regions have a relatively low likelihood, while areas in warmer colors have higher likelihoods.

## HOW IT WORKS

Prior research has shown that, as predicted by gravity models of spatial interaction, the number of co-occurrences between two place names tends to be proportional to the overall frequencies of those place names, and inversely proportional to the distance between them. GRAIL builds on this by computing the likelihood function of an empirically validated model of spatial interaction for a grid of candidate locations in a specified area. 

For every candidate location, the likelihood of the actual pattern of co-occurrences that was observed is computed. These likelihoods are then mapped onto the probabilities of observing a particular place in a particular region, and a heatmap illustrating the likelihoods is drawn to the bin/debug/output/[dataset name] folder. For example, with the default color settings, the odds that the specified toponym's referent can be found in a region that is cyan or warmer is 75%. For each heatmap, a version of the heatmap is generated in which the actual location of the toponym's referent appears as a black square with a green border. Try selecting different datasets to get a sense of the performance of the algorithm under different conditions.

## DEMO DATASETS

Try one of the following datasets by changing the first line of Form1_Load in Form1.cs to one of the following:

DemoSet demoSet = DemoSet.Nunavut;
DemoSet demoSet = DemoSet.Bible;
DemoSet demoSet = DemoSet.USA;

Then, run the application and wait for the computation to complete (will take several minutes). When the computation is complete, a window will pop up instructing you to look in Debug/bin/[dataset name] for the results.
 
### DemoSet.Nunavut

The Legislative Assembly of Nunavut (a Canadian province) publishes the Nunavut Hansard, a verbatim report of legislative proceedings in English and Inuktitut. The present dataset consists of 28 towns and villages, most with populations of less than 1,000, mentioned in the Nunavut Hansard. It is a relatively small corpus (82 MB), but the location of each village can still be estimated with a median error of 366 km. Estimates made using either English or Inuktutit version of the corpus are nearly identical, illustrating the technique's utility in at least one morphologically complex language.

### DemoSet.USA

The Web 1T 5-gram Corpus is a collection of n-grams collected from one terabyte of webpages by Google and released through the Linguistic Data Consortium. The present dataset consists of the 50 largest cities in the United States. Each of these cities is treated as an unknown location, while 725 other cities in the U.S. are treated as known locations. As always, the algorithm is required to estimate each "unknown location" on the basis of co-occurrences with known locations. Performance was quite high on this dataset, showing that high accuracy can be achieved with very large corpora.

### DemoSet.Bible

The ESV Bible Atlas provides latitudes, longitudes, and other geographic information for hundreds of toponyms mentioned in the Bible.  The present dataset consists of 75 towns and villages that are mentioned in at least 4 chapters of the Bible, using the Bible itself as the reference corpus. Compared to most textual corpora used in "big data" applications, the Bible is quite small (approx. 4 megabytes). With such a small amount of data, accuracy is lower than for the Nunavut and United States corpora, but is still much higher than that typically reported in the literature on co-occurrence-based methods for place estimation.