# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import time
import concurrent.futures


def func1():
    while True:
        print("func1")
        time.sleep(3)


def func2():
    while True:
        print("func2")
        time.sleep(3)


# -

executor = concurrent.futures.Prce
f1=executor.submit(func1)
f2=executor.submit(func2)

import os

pid=os.getpid()

# +
import signal, psutil

def kill_child_processes(parent_pid, sig=signal.SIGTERM):
    try:
        parent = psutil.Process(parent_pid)
    except psutil.NoSuchProcess:
        return
    children = parent.children(recursive=True)
    for process in children:
        process.send_signal(sig)


# -

parent = psutil.Process(pid)

parent.children(recursive=True)
