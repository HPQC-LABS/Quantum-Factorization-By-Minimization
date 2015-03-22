# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 12:11:02 2015

@author: Richard
"""

from collections import defaultdict
import itertools
import sqlite3

from semiprime_tools import num_to_factor_num_qubit

DB_NAME = 'semiprimes.db'


INIT_DIMS = [4, 8, 17, 21, 24, 45, 51, 165] + range(10, 520, 10)

PADDING = 3

def spec_tuple_to_id_str((num_fact1, num_fact2, index)):
    ''' Return the convenient string ID '''
    return '{}x{}_{}'.format(num_fact1, num_fact2, str(index).rjust(PADDING, '0'))

def id_str_to_spec_tuple(id_str):
    ''' Return the convenient string ID 

        >>> id_str_to_spec_tuple('60x60_003')    
    '''
    dim_tup, index = id_str.split('_')
    dim1, dim2 = dim_tup.split('x')
    return tuple(map(int, (dim1, dim2, index)))

def iter_rows(tablename):
    ''' Provide an iterable over rows of a table '''
    conn, c = get_connection_and_cursor()
    c.execute('''SELECT * from {}'''.format(tablename))
    rows = c.fetchmany()
    while len(rows):
        for r in rows:
            yield r
        rows = c.fetchmany()

def print_table(tablename, limit=10):
    ''' Print the contents of a table '''
    conn, c = get_connection_and_cursor()
    c.execute('''SELECT * from {} LIMIT {}'''.format(tablename, limit))
    rows = c.fetchall()
    print rows

def get_connection_and_cursor():
    ''' Return a connection to the database '''
    conn = sqlite3.connect(DB_NAME)
    c = conn.cursor()
    return conn, c

def _create_dimension_table():
    ''' Create the table used to store the numbers themselves '''
    conn, c = get_connection_and_cursor()
    try: c.execute('''DROP TABLE dimension''')
    except: pass
    c.execute('''CREATE TABLE dimension (dimension_id INTEGER not null PRIMARY KEY,
                                          num_fact1 int not null,
                                          num_fact2 int not null
                                          )''')
    c.execute('''CREATE UNIQUE INDEX dim_index ON dimension(num_fact1,
                                                           num_fact2);''')
    for dim in INIT_DIMS:
        dim_tuple = (dim, dim)
        c.execute('''INSERT INTO dimension(num_fact1, num_fact2) VALUES (?, ?)''', dim_tuple)
        
    conn.commit()
    conn.close()

def add_new_dimension(dim_tuple):
    ''' Add a new entry to the dimension table '''
    conn, c = get_connection_and_cursor()
    c.execute('''INSERT INTO dimension(num_fact1, num_fact2) VALUES (?, ?)''', dim_tuple)
    conn.commit()

def _create_semiprime_table():
    ''' Create the table used to store the numbers themselves '''
    conn, c = get_connection_and_cursor()
    try: c.execute('''DROP TABLE semiprimes''')
    except: pass
#    c.execute('''CREATE TABLE semiprimes (dimension_id int not null,
#                                          index_ int not null autoincrement,
#                                          value_ text,
#                                          PRIMARY KEY (dimension_id, index_)
#                                          )''')


    c.execute('''CREATE TABLE semiprimes (dimension_id int not null,
                                          index_ int not null,
                                          value_ text not null primary key
                                          )''')

    c.execute('''CREATE UNIQUE INDEX sp_index ON semiprimes(dimension_id,
                                                            index_);''')


    # Save (commit) the changes
    conn.commit()
    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    conn.close()

def dimension_ids_for_semiprimes(semiprimes):
    ''' Return a list of dimension ids, throwing if we have a dimension we
        haven't seen before
    '''
    # Build the map in memory from tuple to id
    tuple_to_id = {}
    for id_, dim1, dim2 in iter_rows('dimension'):
        tuple_to_id[(dim1, dim2)] = id_

    # Now work out the dimension of each semiprime
    dim_tuples = map(num_to_factor_num_qubit, semiprimes)

    return [tuple_to_id[dim_tuple] for dim_tuple in dim_tuples]

def _starting_index_of_dim_id(dim_id):
    ''' Given a dimension id find the next free id '''
    conn, c = get_connection_and_cursor()
    c.execute('''SELECT MAX(semiprimes.index_) WHERE semiprimes.dimension_id IS {}'''.format(dim_id))
    return c.fetchone() + 1

def add_semiprimes_of_same_dimension(semiprimes):
    ''' Given a list of semiprimes with the same dimension, add them to the
        table
    '''
    conn, c = get_connection_and_cursor()
    # fetch the dimension and make sure they're all the same
    dim_ids = dimension_ids_for_semiprimes(semiprimes)
    dim_id = dim_ids.pop()
    for d_id in dim_ids: 
        assert dim_id == d_id
    
    # Now fetch the last index
    c.execute('''SELECT MAX(semiprimes.index_) FROM semiprimes 
                    WHERE semiprimes.dimension_id IS {}'''.format(dim_id))
    start = c.fetchone()[0]
    if start is None:
        start = 1
    else:
        start += 1

    text_semip = map(str, semiprimes)
    c.executemany('''insert into semiprimes(dimension_id, index_, value_) VALUES ({dim_id}, ?, ?)'''.format(dim_id=dim_id),
                                            zip(range(start, start + len(text_semip)), text_semip))
    conn.commit()
    
def add_semiprimes(semiprimes):
    ''' Given an iterable of semiprimes, work out the dimension_id and add them
        to the table
    '''
    dim_ids = dimension_ids_for_semiprimes(semiprimes)
    
    id_to_semiprimes = defaultdict(list)
    for dim_id, sp in itertools.izip(dim_ids, semiprimes):
        id_to_semiprimes[dim_id].append(sp)
    
    for sps in id_to_semiprimes.values():
        add_semiprimes_of_same_dimension(sps)

def semiprimes_to_ids(semiprimes):
    ''' Given a list of semiprimes, find their IDs '''
    conn, c = get_connection_and_cursor()
    
    semi_strs = map(str, semiprimes)
    q_string = ', '.join(['?'] * len(semi_strs))
    c.execute('''SELECT semiprimes.value_, dimension.num_fact1, dimension.num_fact2, semiprimes.index_ FROM semiprimes
                JOIN dimension ON dimension.dimension_id = semiprimes.dimension_id
                WHERE semiprimes.value_ IN (''' + q_string + ')', 
                semi_strs)
    tups = c.fetchall()
    
    # Make a map of semiprime to specification_tuple    
    sp_to_id = {}
    for tup in tups:
        sp_to_id[tup[0]] = tup[1:]

    # Now return the ID strings in the original order
    return [spec_tuple_to_id_str(sp_to_id[sp]) for sp in semi_strs]

def ids_to_semiprimes(str_ids):
    ''' Given a list of IDs, return the corresponding semiprimes '''
    conn, c = get_connection_and_cursor()
    
    spec_tuples = map(id_str_to_spec_tuple, str_ids)
    
#    q_string = ', '.join(['(?, ?, ?)'] * len(spec_tuples))

    c.execute('''SELECT semiprimes.value_, dimension.num_fact1, dimension.num_fact2, semiprimes.index_ FROM semiprimes
                JOIN dimension ON dimension.dimension_id = semiprimes.dimension_id
                ''')
#                WHERE (dimension.num_fact1, dimension.num_fact2, semiprimes.index_) IN (''' + q_string + ')',
#                spec_tuples)

    tups = c.fetchall()
    tups = [(t[0], (t[1], t[2], t[3])) for t in tups]    
    
    # Make a map of semiprime to specification_tuple    
    id_tuple_to_sp = {}
    for tup in tups:
        id_tuple_to_sp[tup[1]] = tup[0]

    # Now return the ID strings in the original order
    return [int(id_tuple_to_sp[id_tup]) for id_tup in spec_tuples]


def _create_factor_table():
    ''' Create the table used to store the numbers themselves '''
    conn, c = get_connection_and_cursor()  
    c.execute('''CREATE TABLE factors
             (dimension_id int, 
             index_ int, 
             factor1 int,
             factor2 int,
             PRIMARY KEY (dimension_id, index_))''')

    # Save (commit) the changes
    conn.commit()
    
    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    conn.close()

def _create_results_table():
    ''' Create the table used to store the numbers themselves '''
    conn, c = get_connection_and_cursor()  
    c.execute('''CREATE TABLE results
             (dimension smallint, 
             index mediumint, 
             qubits_remaining smallint,
             comments text,
             output blob,
             PRIMARY KEY (dimension, index))''')

    # Save (commit) the changes
    conn.commit()
    
    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    conn.close()

if __name__ == '__main__':
#    _create_semiprime_table()
#    _create_dimension_table()
    
    from rsa_constants import RSA100

    
#    dim_tup = num_to_factor_num_qubit(RSA100)
#    add_new_dimension(dim_tup)
#    add_new_dimension((20, 20))
#    add_semiprimes([RSA100, RSA100+1])
#    add_semiprimes([RSA100+3, RSA100+5])
    
#    print_table('dimension')
#    print_table('semiprimes')
    sps = [RSA100+3, RSA100+1]
    ids = semiprimes_to_ids(sps)
    sps2 = ids_to_semiprimes(ids)
    print sps
    print ids
    print sps2
#    print dimension_ids_for_semiprimes([RSA100])