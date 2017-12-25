#!/usr/bin/env python


class verse(Prog_Rsp):

    def __init__(self,Conf,progname,order1={},order2={},silence=False,appendorder=""):
    
        super(HTcount,self).__init__(Conf,progname,order1,order2,silence)
        self.apdord = appendorder
        
    def run(self):
        #for outside use
    
        order = self.gen_order(self.progname,self.data)
        order += self.apdord
        sys.stderr.write(order+"\n")
        return self.run_order(order,self.progname,self.silence)

