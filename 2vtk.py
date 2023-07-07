#!/usr/bin/env python
# encoding: utf-8
'''Convert the binary output of DynEarthSol3D to VTK files.

usage: 2vtk.py [-a -c -m -p -t -h] modelname [start [end [delta]]]]

options:
    -a          save data in ASCII format (default: binary)
    -c          save files in current directory (default: same directory as
                the data files)
    -m          save marker data
    -p          save principal components (s1 and s3) of deviatoric stress
    -t          save all tensor components (default: only 1st/2nd invariants)
    -h,--help   show this help

If 'start' is not provided, start from the 0th frame.
If 'start' is -1, resume previous conversion.
If 'end' is not provided or is -1, end at the last output.
'''

from __future__ import print_function, unicode_literals
import sys, os
import base64, zlib, glob
import numpy as np
from scipy import spatial
from numpy.linalg  import eigh
from Dynearthsol import Dynearthsol

# Save in ASCII or encoded binary.
# Some old VTK programs cannot read binary VTK files.
output_in_binary = True

# Save the resultant vtu files in current directory?
output_in_cwd = False

# Save indivisual components?
output_tensor_components = False

# Save principle stresses
output_principle_stress = False

# Save markers?
output_markers = True

########################
# Is numpy version < 1.8?
eigh_vectorized = True
npversion = np.__version__.split('.')
npmajor = int(npversion[0])
npminor = int(npversion[1])
if npmajor < 1 or (npmajor == 1 and npminor < 8):
    eigh_vectorized = False

class Filter():

  def marker(self,par,x, z, m, t):
    ind = (par.xmin <= x) * (x <= par.xmax) * \
          (par.zmin <= z) * (z <= par.zmax)
    x = x[ind]
    z = z[ind]
    m = m[ind]
    t = t[ind]
    return x, z, m, t

  def node(self,x,z,f):
    ind = (f >= 32 ) * (f <= 34)
    x = x[ind]
    z = z[ind]
    return x, z

########################
# Is numpy version < 1.8?
eigh_vectorized = True
npversion = np.__version__.split('.')
npmajor = int(npversion[0])
npminor = int(npversion[1])
if npmajor < 1 or (npmajor == 1 and npminor < 8):
    eigh_vectorized = False


def main(modelname, start, end, delta):
    prefix = modelname
    if output_in_cwd:
        output_prefix = os.path.basename(modelname)
    else:
        output_prefix = modelname

    des = Dynearthsol(modelname)

    if start == -1:
        vtulist = sorted(glob.glob(modelname + '.*.vtu'))
        lastframe = int(vtulist[-1][(len(modelname)+1):-4]) if vtulist else des.frames[0]
        start = des.frames.index(lastframe) + 1

    if end == -1:
        end = len(des.frames)

    for i in range(start, end, delta):
        frame = des.frames[i]
        nnode = des.nnode_list[i]
        nelem = des.nelem_list[i]
        step = des.steps[i]
        time_in_yr = des.time[i] / (365.2425 * 86400)

        des.read_header(frame)
        suffix = '{0:0=6}'.format(frame)
        print('Converting frame #{0}'.format(suffix), end='\r', file=sys.stderr)

        filename = '{0}.{1}.vtu'.format(output_prefix, suffix)
        fvtu = open(filename, 'w')

        try:
            vtu_header(fvtu, nnode, nelem, time_in_yr, step)

            #
            # node-based field
            #
            fvtu.write('  <PointData>\n')

            # averaged velocity is more stable and is preferred
            try:
                convert_field(des, frame, 'velocity averaged', fvtu)
            except KeyError:
                convert_field(des, frame, 'velocity', fvtu)

            convert_field(des, frame, 'force', fvtu)
            coord0 = des.read_field(frame, 'coord0')
            coord = des.read_field(frame, 'coordinate')
            disp = np.zeros((nnode, 3), dtype=coord.dtype)
            disp[:,0:des.ndims] = coord - coord0
            vtk_dataarray(fvtu, disp, 'total displacement', 3)
            horizon = np.zeros((nnode), dtype=coord.dtype)
            horizon[:] = coord0[:,-1]
            vtk_dataarray(fvtu, horizon, 'horizon', 1)




            '''
            # find nearest neighbour marker of nodes
            markersetname = 'markerset'
            marker_data = des.read_markers(frame, markersetname)
            nmarkers = marker_data['size']
            if nmarkers <= 0:
                raise MarkerSizeError()
            marker_coord = marker_data[markersetname + '.coord']
            marker_mattype = marker_data[markersetname + '.mattype']
            kdtree = spatial.KDTree(marker_coord)
            nn = kdtree.query(coord,1)
            nnmattype = np.zeros((nnode), dtype=marker_mattype.dtype)

            try:
                marker_time = marker_data[markersetname + '.time']
                nnchron = np.zeros((nnode), dtype=marker_time.dtype)
                nnchron[:] = marker_time[nn[1]]
            except:
                pass

            # abjust horizon of sediment node
            horizon_max = max(horizon)
            for j in range(nnode):
                if marker_mattype[nn[1][j]] == 3:
                    horizon[j] = horizon_max
                else:
                    try:
                        nnchron[j] = 0.
                    except:
                        pass


            try:
                vtk_dataarray(fvtu, nnchron, 'chron', 1)
            except:
                pass
            '''

            convert_field(des, frame, 'temperature', fvtu)
            convert_field(des, frame, 'bcflag', fvtu)
            #convert_field(des, frame, 'mass', fvtu)
            #convert_field(des, frame, 'tmass', fvtu)
            #convert_field(des, frame, 'volume_n', fvtu)

            # node number for debugging
            vtk_dataarray(fvtu, np.arange(nnode, dtype=np.int32), 'node number')

            fvtu.write('  </PointData>\n')
            #
            # element-based field
            #
            fvtu.write('  <CellData>\n')

            #convert_field(des, frame, 'volume', fvtu)
            #convert_field(des, frame, 'edvoldt', fvtu)

            convert_field(des, frame, 'mesh quality', fvtu)
            convert_field(des, frame, 'plastic strain', fvtu)
            convert_field(des, frame, 'plastic strain-rate', fvtu)

            strain_rate = des.read_field(frame, 'strain-rate')
            srII = second_invariant(strain_rate)
            vtk_dataarray(fvtu, np.log10(srII+1e-45), 'strain-rate II log10')
            if output_tensor_components:
                for d in range(des.nstr):
                    vtk_dataarray(fvtu, strain_rate[:,d], 'strain-rate ' + des.component_names[d])

            strain = des.read_field(frame, 'strain')
            sI = first_invariant(strain)
            sII = second_invariant(strain)
            vtk_dataarray(fvtu, sI, 'strain I')
            vtk_dataarray(fvtu, sII, 'strain II')
            if output_tensor_components:
                for d in range(des.nstr):
                    vtk_dataarray(fvtu, strain[:,d], 'strain ' + des.component_names[d])

            # averaged stress is more stable and is preferred
            try:
                stress = des.read_field(frame, 'stress averaged')
            except KeyError:
                stress = des.read_field(frame, 'stress')
            tI = first_invariant(stress)
            tII = second_invariant(stress)
            vtk_dataarray(fvtu, tI, 'stress I')
            vtk_dataarray(fvtu, tII, 'stress II')
            if output_tensor_components:
                for d in range(des.ndims):
                    vtk_dataarray(fvtu, stress[:,d] - tI, 'stress ' + des.component_names[d] + ' dev.')
                for d in range(des.ndims, des.nstr):
                    vtk_dataarray(fvtu, stress[:,d], 'stress ' + des.component_names[d])
            if output_principle_stress:
                s1, s3 = compute_principal_stress(stress)
                vtk_dataarray(fvtu, s1, 's1', 3)
                vtk_dataarray(fvtu, s3, 's3', 3)

            convert_field(des, frame, 'density', fvtu)
            convert_field(des, frame, 'material', fvtu)
            convert_field(des, frame, 'viscosity', fvtu)
            effvisc = tII / (srII + 1e-45)
            vtk_dataarray(fvtu, effvisc, 'effective viscosity')

            # element number for debugging
            vtk_dataarray(fvtu, np.arange(nelem, dtype=np.int32), 'elem number')

            # melting mantle
            from scipy import interpolate
            material = des.read_field(frame, 'material')
            temperature = des.read_field(frame, 'temperature')
            connectivity = des.read_field(frame, 'connectivity')
            # Calculate the temperature of element
            ecoord = np.array([coord[connectivity[e,:],:].mean(axis=0) for e in range(nelem)])
            etemp = np.array([temperature[connectivity[e,:]].mean(axis=0) for e in range(nelem)])

            melting = np.zeros(sI.shape)

            # find surface
            bcflag = des.read_field(frame, 'bcflag')
            filter = Filter()
            surfx, surfz = filter.node(coord[:,0],coord[:,1],bcflag)
            orders = np.argsort(surfx)
            surface = np.vstack((surfx[orders],surfz[orders]))
            depth = np.interp(ecoord[:,0], surface[0], surface[1]) - ecoord[:,1]
            pressure = depth * 9.8 * 2900.
            # Hirschmann, 2000  https://doi.org/10.1029/2000GC000070
            # Assumeing solidus is a line between (0 GPa, 1120 C) - (7 GPa, 1800 C)
            #                                                    adibatic  themral gradient 0.3 C/km
            melting[:] = -1000
            ind = material < 2
            melting[ind] = (etemp[ind]-273. + depth[ind]*3.e-4) - (1120 + (680./7.e9)*pressure[ind])
            vtk_dataarray(fvtu, melting, 'melting')

            fvtu.write('  </CellData>\n')

            #
            # node coordinate
            #
            fvtu.write('  <Points>\n')
            convert_field(des, frame, 'coordinate', fvtu)
            fvtu.write('  </Points>\n')

            #
            # element connectivity & types
            #
            fvtu.write('  <Cells>\n')
            convert_field(des, frame, 'connectivity', fvtu)
            vtk_dataarray(fvtu, (des.ndims+1)*np.array(range(1, nelem+1), dtype=np.int32), 'offsets')
            if des.ndims == 2:
                # VTK_ TRIANGLE == 5
                celltype = 5
            else:
                # VTK_ TETRA == 10
                celltype = 10
            vtk_dataarray(fvtu, celltype*np.ones((nelem,), dtype=np.int32), 'types')
            fvtu.write('  </Cells>\n')

            vtu_footer(fvtu)
            fvtu.close()

        except:
            # delete partial vtu file
            fvtu.close()
            os.remove(filename)
            raise

        #
        # Converting marker
        #
        if output_markers:
            # ordinary markerset
            filename = '{0}.{1}.vtp'.format(output_prefix, suffix)
            output_vtp_file(des, frame, filename, 'markerset', time_in_yr, step)

            # hydrous markerset
            if 'hydrous-markerset size' in des.field_pos:
                filename = '{0}.hyd-ms.{1}.vtp'.format(output_prefix, suffix)
                output_vtp_file(des, frame, filename, 'hydrous-markerset', time_in_yr, step)

    print()
    return


def output_vtp_file(des, frame, filename, markersetname, time_in_yr, step):
    fvtp = open(filename, 'w')

    class MarkerSizeError(RuntimeError):
        pass

    try:
        # read data
        marker_data = des.read_markers(frame, markersetname)
        nmarkers = marker_data['size']

        if nmarkers <= 0:
            raise MarkerSizeError()

        # write vtp header
        vtp_header(fvtp, nmarkers, time_in_yr, step)

        # point-based data
        fvtp.write('  <PointData>\n')
        name = markersetname + '.mattype'
        marker_type = marker_data[name]
        vtk_dataarray(fvtp, marker_type, name)
        name = markersetname + '.elem'
        vtk_dataarray(fvtp, marker_data[name], name)
        name = markersetname + '.id'
        vtk_dataarray(fvtp, marker_data[name], name)
        for name in (markersetname + '.time',markersetname + '.z', markersetname + '.distance',markersetname+'.slope'):    
            try:
                vtk_dataarray(fvtp, marker_data[name], name)
            except:
                pass
        fvtp.write('  </PointData>\n')

        # point coordinates
        fvtp.write('  <Points>\n')
        field = marker_data[markersetname + '.coord']
        if des.ndims == 2:
            # VTK requires vector field (velocity, coordinate) has 3 components.
            # Allocating a 3-vector tmp array for VTK data output.
            tmp = np.zeros((nmarkers, 3), dtype=field.dtype)
            tmp[:,:des.ndims] = field
        else:
            tmp = field

        vtk_dataarray(fvtp, tmp, markersetname + '.coord', 3)
        fvtp.write('  </Points>\n')

        vtp_footer(fvtp)
        fvtp.close()

    except MarkerSizeError:
        # delete partial vtp file
        fvtp.close()
        os.remove(filename)
        # skip this frame

    except:
        # delete partial vtp file
        fvtp.close()
        os.remove(filename)
        raise

    return


def convert_field(des, frame, name, fvtu):
    field = des.read_field(frame, name)
    if name in ('coordinate', 'velocity', 'velocity averaged', 'force'):
        if des.ndims == 2:
            # VTK requires vector field (velocity, coordinate) has 3 components.
            # Allocating a 3-vector tmp array for VTK data output.
            i = des.frames.index(frame)
            tmp = np.zeros((des.nnode_list[i], 3), dtype=field.dtype)
            tmp[:,:des.ndims] = field
        else:
            tmp = field

        # Rename 'velocity averaged' to 'velocity'
        if name == 'velocity averaged': name = 'velocity'

        vtk_dataarray(fvtu, tmp, name, 3)
    else:
        vtk_dataarray(fvtu, field, name)
    return


def vtk_dataarray(f, data, data_name=None, data_comps=None):
    if data.dtype in (np.int32,):
        dtype = 'Int32'
    elif data.dtype in (np.single, np.float32):
        dtype = 'Float32'
    elif data.dtype in (np.double, np.float64):
        dtype = 'Float64'
    else:
        raise Error('Unknown data type: ' + name)

    name = ''
    if data_name:
        name = 'Name="{0}"'.format(data_name)

    ncomp = ''
    if data_comps:
        ncomp = 'NumberOfComponents="{0}"'.format(data_comps)

    if output_in_binary:
        fmt = 'binary'
    else:
        fmt = 'ascii'
    header = '<DataArray type="{0}" {1} {2} format="{3}">\n'.format(
        dtype, name, ncomp, fmt)
    f.write(header)
    if output_in_binary:
        header = np.zeros(4, dtype=np.int32)
        header[0] = 1
        a = data.tobytes()
        header[1] = len(a)
        header[2] = len(a)
        b = zlib.compress(a)
        header[3] = len(b)
        f.write(base64.standard_b64encode(header.tobytes()).decode('ascii'))
        f.write(base64.standard_b64encode(b).decode('ascii'))
    else:
        data.tofile(f, sep=' ')
    f.write('\n</DataArray>\n')
    return


def vtu_header(f, nnode, nelem, time, step):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<UnstructuredGrid>
<FieldData>
  <DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">
    {2}
  </DataArray>
  <DataArray type="Float32" Name="CYCLE" NumberOfTuples="1" format="ascii">
    {3}
  </DataArray>
</FieldData>
<Piece NumberOfPoints="{0}" NumberOfCells="{1}">
'''.format(nnode, nelem, time, step))
    return


def vtu_footer(f):
    f.write(
'''</Piece>
</UnstructuredGrid>
</VTKFile>
''')
    return


def vtp_header(f, nmarkers, time, step):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<PolyData>
<FieldData>
  <DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">
    {1}
  </DataArray>
  <DataArray type="Float32" Name="CYCLE" NumberOfTuples="1" format="ascii">
    {2}
  </DataArray>
</FieldData>
<Piece NumberOfPoints="{0}">
'''.format(nmarkers, time, step))
    return


def vtp_footer(f):
    f.write(
'''</Piece>
</PolyData>
</VTKFile>
''')
    return


def first_invariant(t):
    nstr = t.shape[1]
    ndims = 2 if (nstr == 3) else 3
    return np.sum(t[:,:ndims], axis=1) / ndims


def second_invariant(t):
    '''The second invariant of the deviatoric part of a symmetric tensor t,
    where t[:,0:ndims] are the diagonal components;
      and t[:,ndims:] are the off-diagonal components.'''
    nstr = t.shape[1]

    # second invariant: sqrt(0.5 * t_ij**2)
    if nstr == 3:  # 2D
        return np.sqrt(0.25 * (t[:,0] - t[:,1])**2 + t[:,2]**2)
    else:  # 3D
        a = (t[:,0] + t[:,1] + t[:,2]) / 3
        return np.sqrt( 0.5 * ((t[:,0] - a)**2 + (t[:,1] - a)**2 + (t[:,2] - a)**2) +
                        t[:,3]**2 + t[:,4]**2 + t[:,5]**2)


def compute_principal_stress(stress):
    '''The principal stress (s1 and s3) of the deviatoric stress tensor.'''

    nelem = stress.shape[0]
    nstr = stress.shape[1]
    # VTK requires vector field (velocity, coordinate) has 3 components.
    # Allocating a 3-vector tmp array for VTK data output.
    s1 = np.zeros((nelem, 3), dtype=stress.dtype)
    s3 = np.zeros((nelem, 3), dtype=stress.dtype)

    if nstr == 3:  # 2D
        sxx, szz, sxz = stress[:,0], stress[:,1], stress[:,2]
        mag = np.sqrt(0.25*(sxx - szz)**2 + sxz**2)
        theta = 0.5 * np.arctan2(2*sxz, sxx-szz)
        cost = np.cos(theta)
        sint = np.sin(theta)

        s1[:,0] = mag * sint
        s1[:,1] = mag * cost
        s3[:,0] = mag * cost
        s3[:,1] = -mag * sint

    else:  # 3D
        # lower part of symmetric stress tensor
        s = np.zeros((nelem, 3,3), dtype=stress.dtype)
        s[:,0,0] = stress[:,0]
        s[:,1,1] = stress[:,1]
        s[:,2,2] = stress[:,2]
        s[:,1,0] = stress[:,3]
        s[:,2,0] = stress[:,4]
        s[:,2,1] = stress[:,5]

        # eigenvalues and eigenvectors
        if eigh_vectorized:
            # Numpy-1.8 or newer
            w, v = eigh(s)
        else:
            # Numpy-1.7 or older
            w = np.zeros((nelem,3), dtype=stress.dtype)
            v = np.zeros((nelem,3,3), dtype=stress.dtype)
            for e in range(nelem):
                w[e,:], v[e,:,:] = eigh(s[e])

        # isotropic part to be removed
        m = np.sum(w, axis=1) / 3

        p = w.argmin(axis=1)
        t = w.argmax(axis=1)
        #print(w.shape, v.shape, p.shape)

        for e in range(nelem):
            s1[e,:] = (w[e,p[e]] - m[e]) * v[e,:,p[e]]
            s3[e,:] = (w[e,t[e]] - m[e]) * v[e,:,t[e]]

    return s1, s3


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    else:
        for arg in sys.argv[1:]:
            if arg.lower() in ('-h', '--help'):
                print(__doc__)
                sys.exit(0)

    if '-a' in sys.argv:
        output_in_binary = False
    if '-c' in sys.argv:
        output_in_cwd = True
    if '-p' in sys.argv:
        output_principle_stress = True
    if '-t' in sys.argv:
        output_tensor_components = True
    if '-m' in sys.argv:
        output_markers = True

    # delete options
    for i in range(len(sys.argv)):
        if sys.argv[1][0] == '-':
            del sys.argv[1]
        else:
            # the rest of argv cannot be options
            break

    modelname = sys.argv[1]

    if len(sys.argv) < 3:
        start = 0
    else:
        start = int(sys.argv[2])

    if len(sys.argv) < 4 or int(sys.argv[3]) == -1:
        end = -1
    else:
        end = int(sys.argv[3]) + 1

    if len(sys.argv) < 5:
        delta = 1
    else:
        delta = int(sys.argv[4])

    main(modelname, start, end, delta)
