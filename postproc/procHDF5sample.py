#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
this file reads in the hdf5 file that is generated from CAFE and then 
reshapes the dimensions of the data so that it can be read into paraview. it
also creates the corresponding xdmf file
"""

import numpy as np
import os
import h5py

os.chdir('<dirname>')
filbaseT=["<filname1>"] # do not put h5 extension
'''
filbaseT.append("<filname2>") # add as many other files
'''
nFil=len(filbaseT)
for jfil in range(0,nFil):
    filbase = filbaseT[jfil]
    fil=filbase+'.h5'
    f=h5py.File(fil,'r')
    nX=f['Step0/dims'][()]
    nXs=nX[[2,1,0]] # need to swap x and z b/c x is written out first
    dX=f['Step0/VoxelDX'][()]
    dXs=dX[[2,1,0]] # need to swap x and z b/c x is written out first
    tempField=f['Step0/Temperature'][()]
    gID=f['Step0/gID'][()]
    vState=f['Step0/vState'][()]
    IPFz=f['Step0/IPFz'][()]
    IPFx=f['Step0/IPFx'][()]
    IPFy=f['Step0/IPFy'][()]
    f.close()
    gID=np.reshape(gID,nXs)
    tempField=np.reshape(tempField,nXs)
    vState=np.reshape(vState,nXs)
    IPFz=np.reshape(IPFz,np.concatenate((nXs,[3])))
    IPFx=np.reshape(IPFx,np.concatenate((nXs,[3])))
    IPFy=np.reshape(IPFy,np.concatenate((nXs,[3])))
    filh=filbase+'_out.h5'
    fout=h5py.File(filh,'w')
    micro=fout.create_group('micro')
    dsetdim=micro.create_dataset('dims',data=nX)
    dsetdim=micro.create_dataset('VoxelDX',data=dX)
    dsetg=micro.create_dataset('gID',data=gID)
    dsetv=micro.create_dataset('vState',data=vState)
    dsettemp=micro.create_dataset('Temperature',data=tempField)
    dsetipfx=micro.create_dataset('IPFx',data=IPFx)
    dsetipfy=micro.create_dataset('IPFy',data=IPFy)
    dsetipfz=micro.create_dataset('IPFz',data=IPFz)
    fout.close()
    filx=filbase+'_out.xdmf'
    f=open(filx,'w')
    f.write('<?xml version="1.0"?>\n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd"[]>\n')
    f.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
    f.write(' <Domain>\n')
    f.write('  <Grid Name="micro" GridType="Uniform">\n')
    f.write('    <Topology TopologyType="3DCoRectMesh" Dimensions="%d %d %d "></Topology> \n' % (nXs[0],nXs[1],nXs[2]))
    f.write('    <Geometry Type="ORIGIN_DXDYDZ">\n')
    f.write('      <!-- Origin  Z, Y, X -->\n')
    f.write('      <DataItem Format="XML" Dimensions="3">0 0 0</DataItem>\n')
    f.write('      <!-- DxDyDz (Spacing/Resolution) Z, Y, X -->\n')
    f.write('      <DataItem Format="XML" Dimensions="3">%e %e %e</DataItem> \n' % (dXs[0],dXs[1],dXs[2]));
    f.write('    </Geometry>\n')
    
    f.write('    <Attribute Name="IPFx" AttributeType="Vector" Center="Cell">\n')
    f.write('      <DataItem Format="HDF" Dimensions="%d %d %d 3 " NumberType="Float" Precision="4" >\n'\
            % (nXs[0],nXs[1],nXs[2]))
    f.write('%s:micro/IPFx \n'% (filh))
    f.write('      </DataItem>\n')
    f.write('    </Attribute>\n')
    f.write('    <Attribute Name="IPFy" AttributeType="Vector" Center="Cell">\n')
    f.write('      <DataItem Format="HDF" Dimensions="%d %d %d 3 " NumberType="Float" Precision="4" >\n'\
            % (nXs[0],nXs[1],nXs[2]))
    f.write('%s:micro/IPFy \n'% (filh))
    f.write('      </DataItem>\n')
    f.write('    </Attribute>\n')
    f.write('    <Attribute Name="IPFz" AttributeType="Vector" Center="Cell">\n')
    f.write('      <DataItem Format="HDF" Dimensions="%d %d %d 3 " NumberType="Float" Precision="4" >\n'\
            % (nXs[0],nXs[1],nXs[2]))
    f.write('%s:micro/IPFz \n'% (filh))
    f.write('      </DataItem>\n')
    f.write('    </Attribute>\n')
    f.write('    <Attribute Name="gID" AttributeType="Scalar" Center="Cell">\n')
    f.write('      <DataItem Format="HDF" Dimensions="%d %d %d" NumberType="Int" Precision="4" > \n' % (nXs[0],nXs[1],nXs[2]))
    f.write('%s:micro/gID \n' % (filh))
    f.write('      </DataItem>\n')
    f.write('    </Attribute>\n')
    f.write('    <Attribute Name="vState" AttributeType="Scalar" Center="Cell">\n')
    f.write('      <DataItem Format="HDF" Dimensions="%d %d %d" NumberType="Int" Precision="4" > \n' % (nXs[0],nXs[1],nXs[2]))
    f.write('%s:micro/vState \n' % (filh))
    f.write('      </DataItem>\n')
    f.write('    </Attribute>\n')
    f.write('    <Attribute Name="Temperature" AttributeType="Scalar" Center="Cell">\n')
    f.write('      <DataItem Format="HDF" Dimensions="%d %d %d" NumberType="Float" Precision="4" > \n' % (nXs[0],nXs[1],nXs[2]))
    f.write('%s:micro/Temperature \n' % (filh))
    f.write('      </DataItem>\n')
    f.write('    </Attribute>\n')
    f.write('  </Grid>\n')
    f.write(' </Domain>\n')
    f.write('</Xdmf>\n')
    f.close()

#this aggregates all the xdmf files in one such that in paraview one 
# can loop through the data as a time step 
filAgg=filbaseT[0][0:-2]+"_out.xdmf"
fTot=open(filAgg,'w')
fTot.write('<?xml version="1.0"?>\n')
fTot.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd"[]>\n')
fTot.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
fTot.write(' <Domain>\n')
fTot.write('   <Grid Name="CellTime" GridType="Collection" CollectionType="Temporal">\n')
nFil=len(filbaseT)
dt=1.
for jfil in range(0,nFil):
    filbase = filbaseT[jfil]
    fil=filbase+'.h5'
    f=h5py.File(fil,'r')
    filh=filbase+'_out.h5'
    nX=f['Step0/dims'][()]
    nXs=nX[[2,1,0]] # need to swap x and z b/c x is written out first
    dX=f['Step0/VoxelDX'][()]
    dXs=dX[[2,1,0]] # need to swap x and z b/c x is written out first
    f.close()
    tval=jfil*dt
    fTot.write('    <Grid Name="micro" GridType="Uniform">\n')
    fTot.write('     <Time Value="%f" />\n' % (tval))
    fTot.write('      <Topology TopologyType="3DCoRectMesh" Dimensions="%d %d %d "></Topology> \n' % (nXs[0],nXs[1],nXs[2]))
    fTot.write('      <Geometry Type="ORIGIN_DXDYDZ">\n')
    fTot.write('        <!-- Origin  Z, Y, X -->\n')
    fTot.write('        <DataItem Format="XML" Dimensions="3">0 0 0</DataItem>\n')
    fTot.write('        <!-- DxDyDz (Spacing/Resolution) Z, Y, X -->\n')
    fTot.write('        <DataItem Format="XML" Dimensions="3">%e %e %e</DataItem> \n' % (dXs[0],dXs[1],dXs[2]));
    fTot.write('      </Geometry>\n')
    
    fTot.write('      <Attribute Name="IPFx" AttributeType="Vector" Center="Cell">\n')
    fTot.write('        <DataItem Format="HDF" Dimensions="%d %d %d 3 " NumberType="Float" Precision="4" >\n'\
            % (nXs[0],nXs[1],nXs[2]))
    fTot.write('  %s:micro/IPFx \n'% (filh))
    fTot.write('        </DataItem>\n')
    fTot.write('      </Attribute>\n')
    fTot.write('      <Attribute Name="IPFy" AttributeType="Vector" Center="Cell">\n')
    fTot.write('        <DataItem Format="HDF" Dimensions="%d %d %d 3 " NumberType="Float" Precision="4" >\n'\
            % (nXs[0],nXs[1],nXs[2]))
    fTot.write('  %s:micro/IPFy \n'% (filh))
    fTot.write('        </DataItem>\n')
    fTot.write('      </Attribute>\n')
    fTot.write('      <Attribute Name="IPFz" AttributeType="Vector" Center="Cell">\n')
    fTot.write('        <DataItem Format="HDF" Dimensions="%d %d %d 3 " NumberType="Float" Precision="4" >\n'\
            % (nXs[0],nXs[1],nXs[2]))
    fTot.write('  %s:micro/IPFz \n'% (filh))
    fTot.write('        </DataItem>\n')
    fTot.write('      </Attribute>\n')
    fTot.write('      <Attribute Name="gID" AttributeType="Scalar" Center="Cell">\n')
    fTot.write('        <DataItem Format="HDF" Dimensions="%d %d %d" NumberType="Int" Precision="4" > \n' % (nXs[0],nXs[1],nXs[2]))
    fTot.write('  %s:micro/gID \n' % (filh))
    fTot.write('        </DataItem>\n')
    fTot.write('      </Attribute>\n')
    fTot.write('      <Attribute Name="vState" AttributeType="Scalar" Center="Cell">\n')
    fTot.write('        <DataItem Format="HDF" Dimensions="%d %d %d" NumberType="Int" Precision="4" > \n' % (nXs[0],nXs[1],nXs[2]))
    fTot.write('  %s:micro/vState \n' % (filh))
    fTot.write('        </DataItem>\n')
    fTot.write('      </Attribute>\n')
    fTot.write('      <Attribute Name="Temperature" AttributeType="Scalar" Center="Cell">\n')
    fTot.write('        <DataItem Format="HDF" Dimensions="%d %d %d" NumberType="Float" Precision="4" > \n' % (nXs[0],nXs[1],nXs[2]))
    fTot.write('  %s:micro/Temperature \n' % (filh))
    fTot.write('        </DataItem>\n')
    fTot.write('      </Attribute>\n')
    fTot.write('    </Grid>\n')
fTot.write('  </Grid>\n')
fTot.write(' </Domain>\n')
fTot.write('</Xdmf>\n')
fTot.close()




