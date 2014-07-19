#!/usr/bin/python

import math
import random

#-----------------------------------------------------------------------------
def length(n1, n2):
    """Compute the distance between two nodes"""
    return math.sqrt( (n1.x - n2.x)**2 + (n1.y - n2.y)**2 )

#-----------------------------------------------------------------------------
#
class Node( object ):
    def __init__( self, id, x, y ):
        self.id  = id
        self.sid = None
        self.x   = x
        self.y   = y
        self.taken = False
        self.cluster = None

    def __str__( self ):
        out = str( self.id ) + " {} {}".format( self.x, self.y )
        return out

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
    """Compute the total distrance travelled for the given solution"""
    objective = length( solution[-1], solution[0] )
    for index in range(0, len(solution)-1):
        objective += length(solution[index], solution[index+1])
    return objective

#-----------------------------------------------------------------------------
if __name__ == '__main__':
    print("Done")

