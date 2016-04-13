# ImarisReader
Read image and segmentation object data stored in [Imaris](http://www.bitplane.com/) ims files. 

Install
=======

Download the files, unzip, then add the folder to the MATLAB path.

Usage
=====
**Create a reader object to read 'file.ims'.**
    
    fileObj = ImarisReader('file.ims');

**Read image volumes.**
    
    vol = fileObj.DataSet.GetDataVolume(cIdx, tIdx);

* cIdx is the zero-based index for a channel
* tIdx is the zero-based index for a time point

**Read surface vertices and faces.**
    
    vertices = fileObj.Surfaces(1).GetVertices(sIdx);
    faces = fileObj.Surfaces(1).GetTriangles(sIdx);

* sIdx is the zero-based index for a surface

**Read spot positions.**
    
    pos = fileObj.Spots(1).GetPositions;

© 2016, Peter Beemiller (pbeemiller@gmail.com)
