if __name__ == '__main__':

    for i in range(len(outputs)):
        output = ".".join(outputs[i].split('.')[:-1])
        tipsy = nptipsyreader.Tipsy(output)
        tipsy._read_param()
        tipsy._read()
        plots.metalstform(tipsy)                                                                                 
