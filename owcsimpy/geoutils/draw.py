# draw.py

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle
plt.style.use('seaborn-whitegrid')
# We can't expect that the users have latex
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')


import numpy as np
from matplotlib.axes import Axes

def draw(**kwargs):

    from owcsimpy.geoobjects.bases.vector_py import Vector_py as Vector
    from owcsimpy.geoobjects.bases.paramline_py import ParamLine_py as Line
    from owcsimpy.geoobjects.bases.rectplane_py import RectPlane_py as RectPlane
    from owcsimpy.geoobjects.bases.cube_py import Cube_py as Cube
    from owcsimpy.geoobjects.bases.circle_py import Circle_py
    from owcsimpy.geoobjects.models.roomcube_py import RoomCube_py
    from owcsimpy.geoobjects.models.humancube_py import HumanCube_py
    from owcsimpy.geoutils.axes3dmod import Axes3DMod as Axes3D
    from owcsimpy.geoutils.cutils import calcRodriguesMtx as getRodriguesMtx

    from owcsimpy.geoutils.matplotlibmod import pyplotmod as pltmod
    
    assert kwargs is not None, 'Please specify your input'

    listPropHere = ['axes','figure','proj2d',
                    'subplots','nrows','ncols',
                    'vectors','lengths','colors',
                    'pointsources',
                    'baredetectors',
                    'cubes','colorcubevect','enablevect',
                    'lines',
                    'circles','scales',
                    'planes','alphas','facecolors','edgecolors',
                    'models3d','rooms','blockingobjects',
                    'azim','elev',
                    'xlim','ylim','zlim','figsize','aspect']
    
    
    # Separating properties of Line2D and 
    # other items that will be used here
    dictHere = {}
    dictLine = {}
    
    for key,item in kwargs.items():
        if key.lower() in listPropHere:
            dictHere.update({key.lower():item})
        else:
            dictLine.update({key.lower():item})
    
    # Get keys and normalize it to lowercase
    keys = list(dictHere.keys())
    
    # Check one by one
    if 'axes' in keys:
        ax = dictHere.get('axes')

        # will handle 3D only now
        # assert isinstance(ax,Axes3D)
        assert isinstance(ax,Axes)
        
        assert 'figure' in keys

        fig = dictHere.get('figure')

        if 'xlim' in keys:
            ax.set_xlim(dictHere.get('xlim'));
        if 'ylim' in keys:
            ax.set_ylim(dictHere.get('ylim'));
        if 'zlim' in keys:
            ax.set_zlim(dictHere.get('zlim'));

    else:
        figsize = dictHere.get('figsize') if 'figsize' in keys else (5,3)
        aspect = dictHere.get('aspect') if 'aspect' in keys else 1.
        azim = dictHere.get('azim') if 'azim' in keys else -60
        elev = dictHere.get('elev') if 'elev' in keys else 30

        if 'subplots' in keys:
            issubplots = dictHere.get('subplots')
            if issubplots:
                fig = pltmod.figure(figsize=figsize)

                nrows = dictHere.get('nrows') if 'nrows' in keys else 1
                ncols = dictHere.get('ncols') if 'ncols' in keys else 1
                
                axs = []
                for idx in range(nrows*ncols):
                    ax = fig.add_subplot(fig,nrows, ncols, idx+1,  projection='3d')
                    ax.view_init(azim=azim,elev=elev)
                    ax.set_aspect(aspect)
                    ax.set_xlabel(r"$x$");
                    ax.set_ylabel(r"$y$");
                    ax.set_zlabel(r"$z$");
                    axs.append(ax)

                for ax in axs:
                    if 'xlim' in keys:
                        ax.set_xlim(dictHere.get('xlim'));
                    if 'ylim' in keys:
                        ax.set_ylim(dictHere.get('ylim'));
                    if 'zlim' in keys:
                        ax.set_zlim(dictHere.get('zlim'));

            else:
                fig,axs = None,None

        else:
            fig = plt.figure(figsize=figsize)
            ax = Axes3D(fig,aspect=aspect,azim=azim,elev=elev)
        
            if 'xlim' in keys:
                ax.set_xlim(dictHere.get('xlim'));
            if 'ylim' in keys:
                ax.set_ylim(dictHere.get('ylim'));
            if 'zlim' in keys:
                ax.set_zlim(dictHere.get('zlim'));

            ax.set_xlabel(r"$x$");
            ax.set_ylabel(r"$y$");
            ax.set_zlabel(r"$z$");

    
    if 'proj2d' in keys:
        is2D = dictHere.get('proj2d')
        if is2D:
            ax.view_init(azim=-90-1e-3, elev=90-1e-3)
            ax.set_zticks([])
            ax.set_zlabel("")
        
    if 'vectors' in keys or 'pointsources' in keys:
        # Transform into list
        vectors = ([dictHere.get('vectors')] 
             if not isinstance(dictHere.get('vectors'),list) 
             else dictHere.get('vectors'))

        lengths = getElements('lengths',keys,dictHere,1)
        colors = getElements('colors',keys,dictHere,'red')

        # Equalize the length of lists
        if len(lengths) == 1:
            lengths = len(vectors)*lengths

        if len(colors) == 1:
            colors = len(vectors)*colors

        if not(len(lengths)==len(colors)==len(vectors)):
            raise ValueError('Length of elements is not the same.')

        # Draw vector in vector
        for vector,length,color in zip(vectors,lengths,colors):
            if isinstance(vector,Vector):
                dictLine.update({'color':color})
                ax.draw_arrow(vector.refPoint,vector.cartesian,length=length,**dictLine)

    if 'lines' in keys:
        # Transform into list
        lines = ([dictHere.get('lines')] 
             if not isinstance(dictHere.get('lines'),list) 
             else dictHere.get('lines'))

        colors = getElements('colors',keys,dictHere,'black')

        # Equalize the length of lists
        if len(colors) == 1:
            colors = len(lines)*colors

        if not(len(colors)==len(lines)):
            raise ValueError('Length of elements is not the same.')

        # Draw vector in vector
        for line,color in zip(lines,colors):
            if isinstance(line,Line):
                dictLine.update({'color':color})
                # x0,y0,z0 = line.P0
                # x1,y1,z1 = line.P1
                # ax.plot([x0,x1],[y0,y1],[z0,z1],**dictLine)
                ax.draw_arrow(line.P0,line.u,**dictLine)
    
    if 'circles' in keys or 'baredetectors' in keys:
        # Transform into list
        circles = ([dictHere.get('circles')] 
             if not isinstance(dictHere.get('circles'),list) 
             else dictHere.get('circles'))

        # alpha = dictHere.get('alpha') if 'alpha' in keys else 0.4
        # facecolor = dictHere.get('facecolor') if 'facecolor' in keys else 'red'
        # edgecolor = dictHere.get('edgecolor') if 'edgecolor' in keys else 'black'

        alphas = getElements('alphas',keys,dictHere,0.25)
        facecolors = getElements('facecolors',keys,dictHere,'red')
        edgecolors = getElements('edgecolors',keys,dictHere,'black')
        lengths = getElements('lengths',keys,dictHere,1)
        colors = getElements('colors',keys,dictHere,'red')
        scales = getElements('scales',keys,dictHere,1.)

        

        # Equalize the length of lists
        if len(alphas) == 1:
            alphas = len(circles)*alphas

        if len(facecolors) == 1:
            facecolors = len(circles)*facecolors

        if len(edgecolors) == 1:
            edgecolors = len(circles)*edgecolors

        if len(lengths) == 1:
            lengths = len(circles)*lengths

        if len(colors) == 1:
            colors = len(circles)*colors

        if len(scales) == 1:
            scales = len(circles)*scales

        if not(len(alphas)==len(facecolors)==len(edgecolors)==
            len(lengths)==len(colors)==len(circles)):
            raise ValueError('Length of elements is not the same.')

        for circle,alpha,facecolor,edgecolor,length,color,scale in zip(
            circles,alphas,facecolors,edgecolors,lengths,colors,scales):
            if isinstance(circle,Circle_py):
            
                # see https://stackoverflow.com/questions/18228966/how-can-matplotlib-2d-patches-be-transformed-to-3d-with-arbitrary-normals
                pathpatch = Circle((0,0), circle.radius*np.sqrt(scale),
                    facecolor=facecolor,edgecolor=edgecolor,alpha=alpha)
                ax.add_patch(pathpatch)

                path = pathpatch.get_path() #Get the path and the associated transform
                trans = pathpatch.get_patch_transform()

                path = trans.transform_path(path) #Apply the transform

                pathpatch.__class__ = art3d.PathPatch3D #Change the class
                pathpatch._code3d = path.codes #Copy the codes
                pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color    

                verts = path.vertices #Get the vertices in 2D
                pathpatch._segment3d = np.array([np.dot(np.eye(3), (x, y, 0)) + (0, 0, 0) for x, y in verts])

                # Get rotation matrix
                R = getRodriguesMtx(circle.azimuth,np.array([0,0,1]))@getRodriguesMtx(circle.polar,np.array([0,1,0]))

                # Rotate
                for idx,vert in enumerate(pathpatch._segment3d):
                    pathpatch._segment3d[idx] = (R@pathpatch._segment3d[idx].reshape(3,1)).reshape(-1)

                # Translation
                pathpatch._segment3d += circle.ctrPoint
                
                dictLine.update({'color':color})
                ax.draw_arrow(circle.ctrPoint,circle.normalVect,
                    length=length,**dictLine)

    if 'planes' in keys:
        # Transform into list
        planes = ([dictHere.get('planes')] 
             if not isinstance(dictHere.get('planes'),list) 
             else dictHere.get('planes'))

        # alpha = dictHere.get('alpha') if 'alpha' in keys else 0.4
        # facecolor = dictHere.get('facecolor') if 'facecolor' in keys else 'red'
        # edgecolor = dictHere.get('edgecolor') if 'edgecolor' in keys else 'black'

        alphas = getElements('alphas',keys,dictHere,0.4)
        facecolors = getElements('facecolors',keys,dictHere,'red')
        edgecolors = getElements('edgecolors',keys,dictHere,'black')
        lengths = getElements('lengths',keys,dictHere,1)
        colors = getElements('colors',keys,dictHere,'red')

        # Equalize the length of lists
        if len(alphas) == 1:
            alphas = len(planes)*alphas

        if len(facecolors) == 1:
            facecolors = len(planes)*facecolors

        if len(edgecolors) == 1:
            edgecolors = len(planes)*edgecolors

        if len(lengths) == 1:
            lengths = len(planes)*lengths

        if len(colors) == 1:
            colors = len(planes)*colors

        if not(len(alphas)==len(facecolors)==len(edgecolors)==
            len(lengths)==len(colors)==len(planes)):
            raise ValueError('Length of elements is not the same.')

        for plane,alpha,facecolor,edgecolor,length,color in zip(
            planes,alphas,facecolors,edgecolors,lengths,colors):
            assert isinstance(plane,RectPlane)

            if 'enablevect' in keys:
                enablevect = dictHere.get('enablevect')
            else:
                enablevect = True

            # verts = [list(map(lambda vert:
            #     (vert.cartesian[0],vert.cartesian[1],vert.cartesian[2]), 
            #     plane.verts))]
            verts = [list(map(lambda vert: (vert[0],vert[1],vert[2]), 
                plane.verts))]
            col = art3d.Poly3DCollection(verts,alpha=alpha,
                facecolor=facecolor,edgecolor=edgecolor)
            ax.add_collection3d(col)

            if enablevect:
                dictLine.update({'color':color})
                ax.draw_arrow(plane.ctrPoint,plane.normalVect,
                    length=length,**dictLine)

    if 'cubes' in keys:
        # Transform into list
        cubes = ([dictHere.get('cubes')] 
             if not isinstance(dictHere.get('cubes'),list) 
             else dictHere.get('cubes'))

        alphas = getElements('alphas',keys,dictHere,0.4)
        facecolors = getElements('facecolors',keys,dictHere,'red')
        edgecolors = getElements('edgecolors',keys,dictHere,'black')
        lengths = getElements('lengths',keys,dictHere,1)
        colors = getElements('colors',keys,dictHere,'red')

        colorcubevect = dictHere.get('colorcubevect') if 'colorcubevect' in keys else 'red'

        # Equalize the length of lists
        if len(alphas) == 1:
            alphas = len(cubes)*alphas

        if len(facecolors) == 1:
            facecolors = len(cubes)*facecolors

        if len(edgecolors) == 1:
            edgecolors = len(cubes)*edgecolors

        if len(lengths) == 1:
            lengths = len(cubes)*lengths

        if len(colors) == 1:
            colors = len(cubes)*colors

        if not(len(alphas)==len(facecolors)==len(edgecolors)==
            len(lengths)==len(colors)==len(cubes)):
            raise ValueError('Length of elements is not the same.')

        for cube,alpha,facecolor,edgecolor,length,color in zip(
            cubes,alphas,facecolors,edgecolors,lengths,colors):

            assert isinstance(cube,Cube)


            if 'enablevect' in keys:
                enablevect = dictHere.get('enablevect')
            else:
                enablevect = True

            dictLine.update({'color':colorcubevect})
            ax.draw_arrow(cube.ctrPoint,cube.normalVect,
                length=2*length,**dictLine)

            for idx in range(6):
                plane = cube.listPlanes[idx]
                assert isinstance(plane,RectPlane)
                
                verts = [list(map(lambda vert: (vert[0],vert[1],vert[2]), 
                    plane.verts))]
                col = art3d.Poly3DCollection(verts,alpha=alpha,
                    facecolor=facecolor,edgecolor=edgecolor)
                ax.add_collection3d(col)
                if enablevect:
                    dictLine.update({'color':color})
                    ax.draw_arrow(plane.ctrPoint,plane.normalVect,
                        length=length,**dictLine)

    if 'models3d' in keys or 'rooms' in keys or 'blockingobjects' in keys:
        # Transform into list
        models3d = ([dictHere.get('models3d')] 
             if not isinstance(dictHere.get('models3d'),list) 
             else dictHere.get('models3d'))

        # alphas = getElements('alphas',keys,dictHere,0.4)
        facecolors = getElements('facecolors',keys,dictHere,'red')
        edgecolors = getElements('edgecolors',keys,dictHere,'black')
        lengths = getElements('lengths',keys,dictHere,1)
        colors = getElements('colors',keys,dictHere,'red')

        colorcubevect = dictHere.get('colorcubevect') if 'colorcubevect' in keys else 'red'

        # Equalize the length of lists
        # if len(alphas) == 1:
        #     alphas = len(models3d)*alphas

        if len(facecolors) == 1:
            facecolors = len(models3d)*facecolors

        if len(edgecolors) == 1:
            edgecolors = len(models3d)*edgecolors

        if len(lengths) == 1:
            lengths = len(models3d)*lengths

        if len(colors) == 1:
            colors = len(models3d)*colors

        if not(len(facecolors)==len(edgecolors)==
            len(lengths)==len(colors)==len(models3d)):
            raise ValueError('Length of elements is not the same.')

        for model3d,facecolor,edgecolor,length,color in zip(
            models3d,facecolors,edgecolors,lengths,colors):

            assert hasattr(model3d,'listPlanes')

            # if isinstance(model3d,HumanCube_py):
            if not isinstance(model3d,RoomCube_py):
                dictLine.update({'color':color})
                ax.draw_arrow(model3d.ctrPoint,model3d.normalVect,
                    length=2*length,**dictLine)

            # for idx in range(6):
                # plane = model3d.listPlanes[idx]
            for plane in model3d.listPlanes:
                assert isinstance(plane,RectPlane)
                
                verts = [list(map(lambda vert: (vert[0],vert[1],vert[2]), 
                    plane.verts))]
                col = art3d.Poly3DCollection(verts,alpha=0.25*plane.reflectivity,
                    facecolor=facecolor,edgecolor=edgecolor)
                ax.add_collection3d(col)
                if isinstance(model3d,RoomCube_py):
                    dictLine.update({'color':color})
                    ax.draw_arrow(plane.ctrPoint,plane.normalVect,
                        length=length,**dictLine)

    if 'subplots' in keys:
        return fig,axs
    else:
        return fig,ax

def getElements(element,keys,dictHere,defaultVal):
    
    if element in keys:
        listElements = ([dictHere.get(element)] 
            if not isinstance(dictHere.get(element),list) 
            else dictHere.get(element))
    else:
        listElements = [defaultVal]

    return listElements


