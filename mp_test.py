import multiprocessing
from multiprocessing.managers import BaseManager
import stations
from multiprocessing import Process, Manager
import time
import multiprocessing
from multiprocessing.managers import BaseManager
from functools import partial
# class MyManager(BaseManager): pass

def Manager():
    m = BaseManager()
    m.start()
    return m 

class Counter(object):
  def __init__(self):
    self._value = 0

  def update(self, value):
    self._value += value

  def get_value(self):
      return self._value

BaseManager.register('StaLst', stations.StaLst)

# def update(counter_proxy, thread_id):
#   counter_proxy.update(1)
#   print counter_proxy.get_value(), 't%s' % thread_id, \
#     multiprocessing.current_process().name
#   return counter_proxy

# def main():
# def f(counter_p, slst, i):
#     time.sleep(1)
#     print i
#     counter.append(slst[i])
    
def f(counter_p, sta, i):
    # time.sleep(1)
    print i
    counter_p.append(sta)
    
manager = Manager()
counter = manager.StaLst()
SLst = stations.StaLst()
SLst.ReadStaList('station.lst')
pool = multiprocessing.Pool(multiprocessing.cpu_count())
# ADDSLOW = partial(f, datadir=datadir, prefix=prefix, suffix=suffix)
# pool =mp.Pool()
# pool.map(ADDSLOW, self.stations) #make our results with a map call
for i in range(200):
    sta=SLst[i]
    pool.apply_async(func=f, args=(counter, sta, i))
# SLst
pool.close()
pool.join()

# print 'Should be 10 but is %s.' % counter.get_value()