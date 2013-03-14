#!/usr/bin/env python


def flatten_list(list_of_lists):
    """ This is faster than the one-liner version:-
    
    def(flatten): return list(itertools.chain(*list_of_lists)) """

    flat_list = []
    x = flat_list.extend
    for sublist in list_of_lists:
        x(sublist)
    return flat_list
