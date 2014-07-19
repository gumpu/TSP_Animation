#!/usr/bin/python

import math
import random

#-----------------------------------------------------------------------------
def length(n1, n2):
    return math.sqrt( (n1.x - n2.x)**2 + (n1.y - n2.y)**2 )

#-----------------------------------------------------------------------------
class Node( object ):
    def __init__( self, id, x, y ):
        self.id  = id
        self.sid = None
        self.x   = x
        self.y   = y
        self.taken = False
        self.neighbours = None
        self.cluster = None

    def hop( self, nodes ):
        """Hop to the closest neighbour that is not already
        part of the path"""
        min_d = None
        min_n = None
        for n in self.neighbours:
            l = length( self, n )
            if not n.taken:
                if min_d is None:
                    min_d = l
                    min_n = n
                elif l < min_d:
                    min_d = l
                    min_n = n
                else:
                    pass
        if not min_n is None:
            return min_n

        for n in nodes:
            if not n.taken:
                return n

        return None

    def __str__( self ):
        out = str( self.id ) + " {} {}".format( self.x, self.y )
        for i in xrange(5):
            n = self.neighbours[i]
            l = length( n, self )
            out += " {}:{:.5}".format(n.id,l)
        return out

#--------------------------------------------------
def add_neighbours( nodes ):
    number_of_nodes = len( nodes )
    x_axis = [ (n.x, n.id) for n in nodes ]
    y_axis = [ (n.y, n.id) for n in nodes ]

    x_axis.sort()
    y_axis.sort()

    x_index = [0] * number_of_nodes
    y_index = [0] * number_of_nodes

    n = 0
    for i in x_axis:
        x_index[i[1]]=n
        n += 1

    n = 0
    for i in y_axis:
        y_index[i[1]]=n
        n += 1

    for i in xrange(0,number_of_nodes):
        ix = x_index[i]
        iy = y_index[i]

        delta = 10
        go = True
        while go:
            ix_max = ix + delta
            ix_min = ix - delta
            if ix_min < 0:
                ix_max = 2*delta
                ix_min = 0
            elif ix_max >= number_of_nodes:
                ix_min = number_of_nodes-2*delta-1
                ix_max = number_of_nodes-1
            else:
                pass

            iy_max = iy + delta
            iy_min = iy - delta
            if iy_min < 0:
                iy_max = 2*delta
                iy_min = 0
            elif iy_max >= number_of_nodes:
                iy_min = number_of_nodes-2*delta-1
                iy_max = number_of_nodes-1
            else:
                pass

            candidates = [ x_axis[ j ][1] for j in xrange(ix_min, ix_max+1) ]
            neighbours = [ nodes[ y_axis[ j ][1] ] for j in xrange(iy_min, iy_max+1) if y_axis[j][1] in candidates ]

            if len( neighbours ) > 40:
                go = False
                #print i, "OK", len( neighbours )
                nodes[i].neighbours = neighbours[:100]
            else:
                delta = 2*delta
                if delta >= number_of_nodes:
                    delta = number_of_nodes - 1
#-----------------------------------------------------------------------------
def dump( nodes, solution ):
    with open( 'solution.csv', 'w' ) as f: 
        f.write( 'xx,yy,cl\n' )
        for n in solution:
            f.write( str( n.x ) + ',' + str( n.y ) + ',' + str( n.cluster ) + '\n' )
        n = solution[0]
        f.write( str( n.x ) + ',' + str( n.y ) + ',' + str( n.cluster ) + '\n' )

#-----------------------------------------------------------------------------
def total_length( nodes, solution ):
    objective = length( solution[-1], solution[0] )
    for index in range(0, len(solution)-1):
        objective += length(solution[index], solution[index+1])
    return objective

#-----------------------------------------------------------------------------
if __name__ == '__main__':
    random.seed(19671111)

    number_of_nodes = 3000
    nodes = [ Node(i, random.randint(20,80), random.randint(20,80)) for i in xrange(number_of_nodes) ]
    print len( nodes )

    x_axis = [ (n.x, n.id) for n in nodes ]
    y_axis = [ (n.y, n.id) for n in nodes ]

    x_axis.sort()
    y_axis.sort()

    x_index = [0] * number_of_nodes
    y_index = [0] * number_of_nodes

    n = 0
    for i in x_axis:
        x_index[i[1]]=n
        n += 1

    n = 0
    for i in y_axis:
        y_index[i[1]]=n
        n += 1


    if False:
        print x_axis
        print y_axis
        print x_index
        print y_index

    for i in xrange(0,number_of_nodes):
        ix = x_index[i]
        iy = y_index[i]

        delta = 10
        go = True
        while go:
            #print (ix, iy, delta)
            ix_max = ix + delta
            ix_min = ix - delta
            if ix_min < 0:
                ix_max = 2*delta
                ix_min = 0
            elif ix_max >= number_of_nodes:
                ix_min = number_of_nodes-2*delta-1
                ix_max = number_of_nodes-1
            else:
                pass

            iy_max = iy + delta
            iy_min = iy - delta
            if iy_min < 0:
                iy_max = 2*delta
                iy_min = 0
            elif iy_max >= number_of_nodes:
                iy_min = number_of_nodes-2*delta-1
                iy_max = number_of_nodes-1
            else:
                pass

            #print (i, ix_min, ix_max)
            #print (i, ix_min, ix, ix_max)

            candidates = [ x_axis[ j ][1] for j in xrange(ix_min, ix_max+1) ]
            neighbours = [ y_axis[ j ][1] for j in xrange(iy_min, iy_max+1) if y_axis[j][1] in candidates ]
            #print candidates
            if len( neighbours ) > 30:
                go = False
                print i, "OK", len( neighbours )
            else:
                #print i, "====================", delta
                delta = 2*delta
                if delta >= number_of_nodes:
                    delta = number_of_nodes - 1
