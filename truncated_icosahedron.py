import itertools
import numpy as np

def faces(C):
    ''' given list of vertices, find faces, return as list of vertex indices '''

    # for each vertex, find three nearest neighbors
    N = {}
    for ii,c1 in enumerate(C):
        n = []
        for jj,c2 in enumerate(C):
            V = c1 - c2
            D = np.dot(V,V)
            if np.abs(D - 4.0) < .1:
                n.append(jj)
        N[ii] = n
        print ii,len(n),n
    F = []
    centroids = []
    for ii,c1 in enumerate(C):
        for jj in N[ii]:
            for kk in N[jj]:
                if kk!=ii:
                    for ll in N[kk]:
                        D = np.dot(C[ii]-C[ll],C[ii]-C[ll])
                        if ((ll!=jj)&(D<=16.0)):
                            for mm in N[ll]:
                                D = np.dot(C[ii]-C[mm],C[ii]-C[mm])
                                if ((mm!=kk)&(D<=16.0)):
                                    for nn in N[mm]:
                                        D = np.dot(C[ii]-C[nn],C[ii]-C[nn])
                                        if ((nn!=ll)&(D<=16.0)):
                                            if (nn==ii):
                                                Ftmp = np.array([ii,jj,kk,ll,mm])
                                            else:
                                                Ftmp = np.array([ii,jj,kk,ll,mm,nn])
                                            V = 0.0
                                            n = 0.0
                                            for rr in Ftmp:
                                                V += C[rr]
                                                n += 1.0
                                            CC = V / n
                                            good = True
                                            for CCC in centroids:
                                                if (np.dot(CCC-CC,CCC-CC) < 1.0):
                                                    good = False
                                            if (good):
                                                F.append(Ftmp)
                                                centroids.append(CC)

    return F, np.array(centroids)


def edges(C):
    ''' take list of verticies, find edges (pairs of v spearated by exactly 2) '''

    E = []
    for ii,c1 in enumerate(C):
        for jj,c2 in enumerate(C[ii+1:]):
            V = c1 - c2
            D = np.dot(V,V)
            if D==4.0:
                E.append([c1,c2])

    return E

def makevertices():
    '''generate a list of coordinates for the verticies of the truncated icosahedron'''
    phi = 0.5*(1+np.sqrt(5))
    c1 = (0,1,3*phi)
    c2 = (2,(1+2*phi),phi)
    c3 = (1,2+phi,2*phi)
    combos1 = list(itertools.product((1,-1),repeat=2))
    for i in range(len(combos1)):
        combos1[i] = (1,)+combos1[i]
    combos23 = list(itertools.product((1,-1),repeat=3))
    coords = []
    for i in combos1:
        coords.append(np.matrix(map(lambda x,y:x*y,c1,i)).transpose()) #column vectors
    for i in combos23:
        coords.append(np.matrix(map(lambda x,y:x*y,c2,i)).transpose())
        coords.append(np.matrix(map(lambda x,y:x*y,c3,i)).transpose())
    #permutation matrices
    P1 = np.matrix([[0,0,1],[1,0,0],[0,1,0]])
    P2 = np.matrix([[0,1,0],[0,0,1],[1,0,0]])
    for i in coords[:]:
        coords.append((P1*i))
        coords.append((P2*i))
        #coords = [tuple(i.transpose().tolist()[0]) for i in coords]
    coords = np.array([np.array(i.transpose().tolist()[0]) for i in coords])
    return coords


class TruncatedIcosahedron(object):

    def __init__(self):
        self.v = makevertices()
        self.faces, self.centroids = faces(self.v)
        self.phi = 0.5*(1+np.sqrt(5))
        self.r = np.sqrt(9 * self.phi + 10)
        self.rotate_pentagon_to_northpole()
        top_v = self.faces[2][np.argmax(self.v[self.faces[2], 2])]
        self.ref_azimuth = 0
        self.ref_azimuth = self.get_azimuth(top_v)

    def find_adjacent(self, face_index, edge):
        adj = np.where(map(set(edge).issubset, self.faces))[0]
        return adj[adj != face_index][0]

    def get_azimuth(self, v_index):
        az = - self.ref_azimuth + np.arctan2(self.v[v_index][1], self.v[v_index][0]) + np.pi
        if az < 0:
            az += 2*np.pi
        print "Vertex %d, Az %.2f" % (v_index, np.degrees(az))
        return az

    def find_face_south(self, face_index):
        face = self.faces[face_index]
        edge_south = face[np.argsort(self.v[face, 2])[:2]]
        print "Edge south", edge_south
        return self.find_adjacent(face_index, edge_south)

    def find_face_north(self, face_index):
        face = self.faces[face_index]
        edge_south = face[np.argsort(self.v[face, 2])[-2:]]
        print "Edge north", edge_south
        return self.find_adjacent(face_index, edge_south)

    def find_face_se(self, face_index):
        face = self.faces[face_index]
        edge_se = [np.argmax(map(self.get_azimuth, face))]
        nextc = edge_se[0]+1
        if nextc >= len(face):
            nextc = 0
        edge_se.append(nextc)
        if self.v[face[nextc],2] > self.v[face[edge_se[0]-1], 2]:
            edge_se[-1] = edge_se[0]-1
        edge_se = [face[ed] for ed in edge_se]
        print "Edge SE", edge_se
        return self.find_adjacent(face_index, edge_se)

    def find_face_ne(self, face_index):
        face = self.faces[face_index]
        edge_se = [np.argmax(map(self.get_azimuth, face))]
        nextc = edge_se[0]+1
        if nextc > len(face):
            nextc = 0
        edge_se.append(nextc)
        if self.v[face[nextc],2] < self.v[face[edge_se[0]-1], 2]:
            edge_se[-1] = edge_se[0]-1
        edge_se = [face[ed] for ed in edge_se]
        print "Edge NE", edge_se
        return self.find_adjacent(face_index, edge_se)

    def rotate_pentagon_to_northpole(self):
        rot_angle = np.arccos(np.dot(self.centroids[2]/np.linalg.norm(self.centroids[2]), np.array([0,0,1])))

        rotmatrix = np.array([[ 1,       0          ,          0        ],
                             [ 0, np.cos(rot_angle), -np.sin(rot_angle)],
                             [ 0, np.sin(rot_angle),  np.cos(rot_angle)]])

        self.v = np.array([np.dot(rotmatrix, Crow) for Crow in self.v])
        self.faces, self.centroids = faces(self.v)

        for x in self.centroids[2][0:2]:
            assert np.abs(x) < 1e-10
