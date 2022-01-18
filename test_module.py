import sys 

def read_args():
    '''
    Read arguments and corresponding values storing them in a 
    dictionary

    DESCRIPTION
    -----------
    Arguments should be specified in the form arg=value.
    If no value is specified, the value of that argument will be True.
    Even if no argument is specified, there will always an argument 
    called "name" containing the name of the script itself

    RETURNS
    -------
    arg_dict: dictionary
        Dictionary containing argument:value

    HISTORY
    -------
    2021 10 13, Stefano Rapisarda (Uppsala), cretion date
    '''
    args=sys.argv

    arg_dict = {}
    for i,arg in enumerate(args):
        if i == 0:
            div = ['name',arg]
        else:
            if ':' in arg:
                div = arg.split(':')
            elif '=' in arg:
                div = arg.split('=')
            else:
                div = [arg.strip(),True]
        if type(div[0]) == str: div[0]=div[0].strip()
        if type(div[1]) == str: div[1]=div[1].strip()
        arg_dict[div[0]] = div[1] 

    return arg_dict