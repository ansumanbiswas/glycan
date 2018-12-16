''' code for representing the glycan structure, given the input nodes (i.e. the individual sugars)'''
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
    def __init__(self, rootName):
        self.__node = Node(rootName)
        self.__leftChild = None
        self.__rightChild = None

    def add_left_child(self, nodeName):
        if self.__leftChild is None:
            self.__leftChild = TreeNode(nodeName)

    def add_right_child(self, nodeName):
        if self.__rightChild is None:
            self.__rightChild = TreeNode(nodeName)

    def get_left_child(self):
        return self.__leftChild

    def get_right_child(self):
        return self.__rightChild

    def has_left_child(self):
        return self.__leftChild is not None

    def has_right_child(self):
        return self.__rightChild is not None

    def __str__(self):
        repStr = ''
        if self.__leftChild is not None:
            repStr = '(%s)' % self.__leftChild
        repStr = repStr + '%s' % self.__node
        if self.__rightChild is not None:
            repStr = repStr + '(%s)' % self.__rightChild
        return repStr

    def search(self, nodeName):
        res = []
        if self.__node == Node(nodeName):
            res.append(self)
        if self.__leftChild is not None:
            res = res + self.__leftChild.search(nodeName)
        if self.__rightChild is not None:
            res = res + self.__rightChild.search(nodeName)
        return res


def attach_new_node( root, basenode, newnode, attach_to='any', strategy='random'):
    assert isinstance(root, TreeNode)
    assert strategy in ['random', 'greedy']
    assert attach_to in ['left', 'right', 'any']
    hits = root.search(basenode)

    if len(hits) > 0:
        filtered = []
        for h in hits:
            if attach_to == 'left' and not h.has_left_child():
                filtered.append(h)
            elif attach_to == 'right' and not h.has_right_child():
                filtered.append(h)
            elif attach_to == 'any' and (not h.has_right_child() or not h.has_right_child()):
                filtered.append(h)
        if len(filtered) > 0:
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
     print( len(freq) )
'''

