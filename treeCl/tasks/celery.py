#!/usr/bin/env python
from __future__ import absolute_import

import os
from celery import Celery


try:
    homedir = os.getenv('HOME')
    with open(os.path.join(homedir, '.broker', 'connection')) as fl:
        broker_conn = fl.readline().strip()
    with open(os.path.join(homedir, '.backend', 'connection')) as fl:
        backend_conn = fl.readline().strip()

    app = Celery('treeCl.tasks',
                 broker='redis://{}'.format(broker_conn),
                 backend='redis://{}'.format(backend_conn),
                 include=['treeCl.tasks.tasks'])

    # Optional configuration, see the application user guide.
    app.conf.update(
        CELERY_TASK_RESULT_EXPIRES=36000,
    )

except IOError:
    app = Celery('treeCl.tasks',
                 include=['treeCl.tasks.tasks'])


if __name__ == '__main__':
    app.start()
