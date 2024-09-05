// StackProfileData
// This ImageJ macro gets the profile of all slices in a stack
// and writes the data to the Results table, one column per slice.

macro "Stack profile Data" {
     if (!(selectionType()==0 || selectionType==5 || selectionType==6))
       exit("Line or Rectangle Selection Required");
     setBatchMode(true);

     run("Plot Profile");
     Plot.getValues(x, y);
     run("Clear Results");
     for (i=0; i<x.length; i++)
         setResult("x", i, x[i]);
     close();

     n = nSlices;
     for (slice=1; slice<=n; slice++) {
         showProgress(slice, n);
         setSlice(slice);
         profile = getProfile();
         sliceLabel = toString(slice);
         sliceData = split(getMetadata("Label"),"\n");
         if (sliceData.length>0) {
             line0 = sliceData[0];
             if (lengthOf(sliceLabel) > 0)
                 sliceLabel = sliceLabel+ " ("+ line0 + ")";
         }
         for (i=0; i<profile.length; i++)
             setResult(sliceLabel, i, profile[i]);
     }
     setBatchMode(false);
     updateResults;
}
