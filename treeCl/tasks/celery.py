#!/usr/bin/env python
from __future__ import absolute_import

import os

from celery import Celery


homedir = os.getenv('HOME')

try:
    with open(os.path.join(homedir, '.broker', 'connection')) as fl:
        broker_conn = fl.readline().strip()
    with open(os.path.join(homedir, '.backend', 'connection')) as fl:
        backend_conn = fl.readline().strip()
except IOError, e:
    import sys

    sys.stderr.write('FATAL ERROR: no redis connection available - exiting\n')
    sys.stderr.flush()
    sys.exit()

app = Celery('treeCl.tasks',
             broker='redis://{}'.format(broker_conn),
             backend='redis://{}'.format(backend_conn),
             include=['treeCl.tasks.tasks'])

# Optional configuration, see the application user guide.
app.conf.update(
    CELERY_TASK_RESULT_EXPIRES=36000,
)

if __name__ == '__main__':
    app.start()
