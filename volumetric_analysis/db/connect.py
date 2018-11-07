"""
connect.py


"""
import MySQLdb

#Brittin modules

def default(db):
    return  MySQLdb.connect(read_default_file="~/.my.cnf",     
                            db=db)


def connect(host,user,passwd,db):
    return MySQLdb.connect(host=host, 
                           port=3306, 
                           user=user, 
                           passwd=passwd, 
                           db=db)
    
