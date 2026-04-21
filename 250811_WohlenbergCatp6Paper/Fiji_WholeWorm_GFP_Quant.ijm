
// This macro batch processes all the files in a folder and any
// subfolders in that folder. In this example, it runs... 
//For other kinds of processing,
// edit the processFile() function at the end of this macro.
// Probably need to edit the "open" path in line 9 w/ s/t on your machine that contains TIFFs.
// 		Location doesn't matter b/c you will select the actual location after starting the run.

	"BatchProcessFolders"
   requires("1.33s"); 
   dir = getDirectory("Choose a Directory ");
   setBatchMode(true);
   count = 0;
   countFiles(dir);
   n = 0;
   processFiles(dir);
   //print(count+" files processed");
   
   function countFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/"))
              countFiles(""+dir+list[i]);
          else
              count++;
      }
  }

   function processFiles(dir) {
      list = getFileList(dir);
      for (i=0; i<list.length; i++) {
          if (endsWith(list[i], "/"))
              processFiles(""+dir+list[i]);
          else {
             showProgress(n++, count);
             path = dir+list[i];
             processFile(path);
          }
      }
  }

  function processFile(path) {
  	//change file type if not working with .tif
       if (endsWith(path, ".tif")) {
           open(path);
         
function QuantifyParticles(){
run("ROI Manager...");
run("Set Scale...", "distance=1.2403 known=1 pixel=1 unit=micron");
G =getTitle();
selectWindow(G);
print (G);
run("Duplicate...", " ");
D =getTitle();
selectWindow(D);
//print (D);
run("Auto Threshold", "method=Li white");
run("Set Measurements...", "area mean add redirect=["+G+"] decimal=0");
run("Analyze Particles...", "size=5000-Infinity show=Outlines display summarize add");
selectWindow(G);
roiManager("Select", 0);
run("Flatten");
roiManager("Delete");
T= "flat";
saveAs("Tiff");
run("Close All");
}
QuantifyParticles();
}
}
