#!/usr/bin/env python
from __future__ import absolute_import

import os
from celery import Celery


class MockInspect(object):
    def active(self):
        return None

class MockControl(object):
    def inspect(self):
        return MockInspect()

class MockApp(object):
    control = MockControl()

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
    app = MockApp()


if __name__ == '__main__':
    if isinstance(app, MockApp):
        import sys
        sys.stderr.write('ERROR: no redis connection available\n')
        sys.stderr.flush()
        sys.exit()
    app.start()
