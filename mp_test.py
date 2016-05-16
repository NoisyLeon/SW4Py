from multiprocessing import Process, Manager, Lock
import os
import time
from multiprocessing import Process, Value, Array
def f(d, lst, l):
    time.sleep(1)
    d[1] = '1'
    d['2'] = 2
    d[0.25] = None
    l.acquire()
    lst.reverse()
    print lst, os.getpid()
    l.release()
l=Lock()
# def f(n, a, l):
#     n.value = 3.1415927
#     time.sleep(1)
#     l.acquire()
#     for i in range(len(a)):
#         a[i] = -a[i]
#     print arr[:], os.getpid()
#     l.release()
    # print num.value
    
   
manager = Manager()

d = manager.dict()
lst = manager.list(range(10))
num = Value('d', 0.0)
arr = Array('i', range(10))
p=[]
for i in range(5):
    p.append( Process(target=f, args=(d, lst, l)))
# p2=Process(target=f, args=(d, lst, l))
# p3=Process(target=f, args=(d, lst, l))
for i in range(5):
    p[i].start()
for i in range(5):
    p[i].join()

# p3.start()
# p1.start()
# p2.start()
# p1.join()
# p2.join()
# p3.join()
    # p[i].start
# for i in range(5):
    # p[i].join()

    # print d
    # print l