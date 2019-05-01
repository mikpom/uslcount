from pkg_resources import resource_filename as pkg_file
import json
import numpy as np

pred_feats = {'monoexonic' : ['unq_frac', 'nearbyRPK'],
              'multiexonic': ['exo2intr', 'jctF', 'unq_frac', 'nearbyRPK'],
              'multicovered' : ['jctF', 'unq_frac', 'nearbyRPK']}

def get_leafs(tree):
    leafs = []
    def recur(node):
        for descendant in node['left'], node['right']:
            if not descendant is None:
                # could check only one of them; #justToMakeSure
                if descendant['left'] is None and descendant['right'] is None:
                    leafs.append(descendant)
                recur(descendant)
    recur(tree)
    return leafs

# Load decision trees
dtrees = {}
probs = []
for gncls in pred_feats:
    with open(pkg_file('uslcount', 'models/'+gncls+'.json'), 'rt') as js:
        tree = json.load(js)
        dtrees[gncls] = tree
        for leaf in get_leafs(tree):
            probs.append(leaf['value'][0]/leaf['n_node_samples'])
minprob = np.min(probs)
maxprob = np.max(probs)

# Function idea taken from this Stack Overflow answer
# https://stackoverflow.com/a/50717952/6762775
def export_dict(tree, max_depth=None) :
    """
    Export a decision tree in dict format.
    Parameters
    ----------
    decision_tree : decision tree classifier
        The decision tree to be exported
    max_depth : int, optional (default=None)
        The maximum depth of the representation. If None, the tree is fully
        generated.
    Returns
    -------
    a dictionary of the format <tree> := {
        'node_id': <int>,
        'isleaf' : <bool>,
        'feature': <int> | <string>,
        'threshold' : <float>,
        'impurity' : <float>,
        'n_node_samples' : <int>,
        'left' : <tree>,
        'right' : <tree>,
        'value' : [<int>],
    }

    For leaf nodes these fields are abscent: feature, threshold, impurity
    """
    tree_ = tree.tree_

    # i is the element in the tree_ to create a dict for
    def recur(i, depth=0) :
        i = int(i)
        if max_depth is not None and depth > max_depth :
            return None

        # Common attributes
        nd = {'node_id':i,
              'value':tree_.value[i][0].tolist(),
              'n_node_samples' : int(tree_.n_node_samples[i])}

        if tree_.children_left[i] == -1 : # leaf
            nd.update({'left':None,
                       'right':None,
                       'isleaf':True})
        else: # non-leaf
            nd.update({'feature':int(tree_.feature[i]),
                       'threshold':float(tree_.threshold[i]),
                       'impurity' : float(tree_.impurity[i]),
                       'left'  : recur(tree_.children_left[i],  depth + 1),
                       'right' : recur(tree_.children_right[i], depth + 1),
                       'isleaf':False})
            if i == 0: # root
                nd.update({'classes':tree.classes_.tolist()})
        return nd
            
    return recur(0)

def predict_proba(X, dtree):
    """
    Outputs list of probabilities of beloning to classes for a given array of
    features and a decision tree.

    X - array of values. Length should be equal to 
        number of features in decision tree. 
    dtree - decision tree in a form of a recursive dictionary.

    """
    cur_node = dtree
    
    # traverse untill leaf
    while not cur_node['isleaf']:
        if X[cur_node['feature']]<=cur_node['threshold']:
            cur_node = cur_node['left']
        else:
            cur_node = cur_node['right']

    # probabilites
    _pr = [v/cur_node['n_node_samples'] for v in cur_node['value']]

    # make a dictionary
    pr = dict(zip(dtree['classes'], _pr))
    return pr
