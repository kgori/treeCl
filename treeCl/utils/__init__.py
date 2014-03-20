# import fileIO
def flatten_list(list_):
    newlist = list()
    x = newlist.extend
    for sublist in list_:
        x(sublist)
    return newlist

def regex_search_extract(search_attempt):
    return (search_attempt.group() if search_attempt else None)
