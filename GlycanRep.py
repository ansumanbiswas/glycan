''' 
code for representing the glycan structure, given the input nodes (i.e. the individual sugars). 
Each node can have at most 2 branches (left child & right child)
Define class to describe node
Define class to describe: adding child node, get identity of child node, check if node has children.
Define function for attachment strategy
Define rule set
Define simulation protocol
'''

import pandas as pd
import numpy as np

__all__ = ['Node', 'TreeNode', 'attach_new_node']  # each treeNode can either have a node (i.e, child) or be null (having no child)


class Node:      ## WHAT IS THIS CLASS DOING ?
    def __init__(self, id):
        self.__name = id

    def __str__(self):
        return '%s' % self.__name

    def __eq__(self, rhs):
        if isinstance(rhs, Node):
            return self.__name == rhs.__name

    def name(self):
        return self.__name


class TreeNode:
    def __init__(self, rootName):        # defining the name of each node
        self.__node = Node(rootName)
        self.__leftChild = None
        self.__rightChild = None

    def add_left_child(self, nodeName):  # adding left child node
        if self.__leftChild is None:
            self.__leftChild = TreeNode(nodeName)

    def add_right_child(self, nodeName):  # adding right child node
        if self.__rightChild is None:
            self.__rightChild = TreeNode(nodeName)

    def get_left_child(self):          # getting the identity of left child node
        return self.__leftChild

    def get_right_child(self):         # getting the identity of right child node
        return self.__rightChild

    def has_left_child(self):          # check if the selected node has left child node
        return self.__leftChild is not None

    def has_right_child(self):         # check if the selected node has right child node
        return self.__rightChild is not None

    def __str__(self):                # check if the selected node has left & right children; if yes, return the child node(s)
        repStr = ''
        if self.__leftChild is not None:
            repStr = '(%s)' % self.__leftChild
        repStr = repStr + '%s' % self.__node
        if self.__rightChild is not None:
            repStr = repStr + '(%s)' % self.__rightChild
        return repStr

    def search(self, nodeName):       # if a specified node is present in the tree, return the nodes & its children
        res = []
        if self.__node == Node(nodeName):
            res.append(self)
        if self.__leftChild is not None:
            res = res + self.__leftChild.search(nodeName)
        if self.__rightChild is not None:
            res = res + self.__rightChild.search(nodeName)
        return res


def attach_new_node( root, basenode, newnode, attach_to='any', strategy='random'):  # strategy of attaching child nodes 'newnode' to parent node 'basenode'  
    assert isinstance(root, TreeNode)      # this checks that the specified root & intermediate nodes are present in the tree
    assert strategy in ['random', 'greedy']  # strategy of attachment; random: attach to a random subset of the specified nodes;
    # greedy: attach to any one of the multiple occurrences of a specified node 
    assert attach_to in ['left', 'right', 'any']
    hits = root.search(basenode)  # check if the basenode is already present in tree, for further attachment to proceed

    if len(hits) > 0:   ## this section attaches childnode to parent node if empty slots are available
        filtered = []
        for h in hits:
            if attach_to == 'left' and not h.has_left_child():
                filtered.append(h)
            elif attach_to == 'right' and not h.has_right_child():
                filtered.append(h)
            elif attach_to == 'any' and (not h.has_right_child() or not h.has_right_child()):
                filtered.append(h)
        if len(filtered) > 0:  ## stop after 1st attachment for greedy; else try further attachment
            if strategy == 'greedy':
                idx = 0
            else:
                idx = np.random.choice(list(range(len(filtered))))
            if attach_to == 'left':
                filtered[idx].add_left_child(newnode)
                return True
            elif attach_to == 'right':
                filtered[idx].add_right_child(newnode)
                return True
            elif np.random.randint(1,100) % 2 == 1:
                filtered[idx].add_left_child(newnode)
                return True
            else:
                filtered[idx].add_right_child(newnode)
                return True
    return False



'''
if __name__ == "__main__" :
     rules = [ ("N1", "N2", "left"),
               ("N2", "N5", "left"),
               ("N1", "N3", "right"),
               ("N3", "N4", "left"),
               ("N5", "N3", "left"),
               ("N5", "N1", "right") ]

     # simulation
     N = 10000
     freq = {}
     for i in range(N):
         root = TreeNode("N1")
         np.random.shuffle(rules)
         for b, a, r in rules:
             attach_new_node(root, b, a, attach_to=r)
         tree = "%s" % root
         if tree not in freq:
             freq[tree] = 0
         freq[tree] = freq[tree] + 1
     print( len(freq) )     ## check the freq of occurrence of the different trees
'''

