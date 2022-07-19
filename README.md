# colocAnalysis
Microglia-SST colocalization code for Gesuita et al., 2022, Cell Reports.

%%%%%
To automatically count microglia putative interactions with SST+ processes, we first segmented SST+ cell bodies and SST+ processes separately. While cell bodies are blob-like, processes are tubular-like. We enhanced blob-like and tubular-like signals using LoG (Kong et al., 2013) and Jerman vesselness (Jerman et al., 2015) methods respectively. Then, we segmented processes with Otsu thresholding (Otsu, 1979). To improve cell body detection, we then performed connected component analysis of all segmented pieces and extracted structural features (such as roundness, area, mean intensity, etc.)  and clustered into two classes using k-means method (Hartigan and Wong, 1979). We manually determined the class boundary to separate cell bodies from other false-positives and trained a support vector machine (Scholkopf and Smola, 2018) classifier. To segment microglia processes, we used adaptive k-means clustering (k-means and Otsu gives similar results under Euclidean distance metric) (Bishop, 2006). After the segmentation of all structures separately from maximum intensity projected images, we turned these 2-dimensional segmentations into 3-dimensions by thresholding each z-profile that falls under a positively segmented pixel. The final colocalization were performed by simple logical operators in 3-dimensions afterwards. Analysis was performed in MATLAB (MathWorks).
%%%%%

Download Bioformats Matlab Library and Put in this same folder with this code (https://downloads.openmicroscopy.org/bio-formats/5.3.4/artifacts/bfmatlab.zip)

Run colocAnalysis.m
