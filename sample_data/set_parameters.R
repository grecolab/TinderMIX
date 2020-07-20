# parameters for phenodata
dose_index = 2
time_point_index = 3

# #parameter for anova
# pvalAnova = 0.01
# pvalAnova.adj.method = "none"

#parameter for fitting
gridSize = 50
pvalFitting = 0.01
pvalFitting.adj.method = "fdr"
modelSelection = 1:3
logScale = TRUE

# parameter BMD identification
activity_threshold = 0.1
BMD_response_threshold = 0.5
mode = "most_left"
timeLabels =  c("Late","MiddleL","MiddleE","Early")
doseLabels = c("Sensitive","Intermediate","Resilient")
nTimeInt = 4
nDoseInt = 3

#parameter clusterings
nClust = c(5,10,15,20,25)
method="pearson"
hls.method = "ward"

#parameter enrichment
corrType = "fdr"
type_enrich="KEGG"
org_enrich = "rnorvegicus"
pth = 0.05
sig = FALSE
mis = 0
only_annotated=FALSE
