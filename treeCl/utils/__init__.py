# import fileIO
def flatten_list(list_):
    newlist = list()
    x = newlist.extend
    for sublist in list_:
        x(sublist)
    return newlist
