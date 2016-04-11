"""
http://stackoverflow.com/a/28708604

MySQLdb session module for CherryPy by Ken Kinder <http://kenkinder.com/>

Version 0.3, Released June 24, 2000.

Copyright (c) 2008-2009, Ken Kinder
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

    * Neither the name of the Ken Kinder nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import MySQLdb
import cPickle as pickle
import cherrypy
import logging
import threading

__version__ = '0.2'

logger = logging.getLogger('Session')

class MySQLSession(cherrypy.lib.sessions.Session):
    ##
    ## These can be over-ridden by config file
    table_name = 'web_session'
    connect_arguments = {}

    SCHEMA = """create table if not exists %s (
            id varchar(40),
            data text,
            expiration_time timestamp
        ) ENGINE=InnoDB;"""

    _database = None

    def __init__(self, id=None, **kwargs):
        logger.debug('Initializing MySQLSession with %r' % kwargs)
        for k, v in kwargs.items():
            setattr(MySQLSession, k, v)

        self.db = self.get_db()
        self.cursor = self.db.cursor()

        super(MySQLSession, self).__init__(id, **kwargs)

    @classmethod
    def get_db(cls):
        ##
        ## Use thread-local connections
        local = threading.local()
        if hasattr(local, 'db'):
            return local.db
        else:
            logger.debug("Connecting to %r" % cls.connect_arguments)
            db = MySQLdb.connect(**cls.connect_arguments)
            cursor = db.cursor()
            cursor.execute(cls.SCHEMA % cls.table_name)
            db.commit()
            local.db = db

            return db

    def _load(self):
        logger.debug('_load %r' % self)
        # Select session data from table
        self.cursor.execute('select data, expiration_time from %s '
                            'where id = %%s' % MySQLSession.table_name, (self.id,))
        row = self.cursor.fetchone()
        if row:
            (pickled_data, expiration_time) = row
            data = pickle.loads(pickled_data)

            return data, expiration_time
        else:
            return None

    def _save(self, expiration_time):
        logger.debug('_save %r' % self)
        pickled_data = pickle.dumps(self._data)

        self.cursor.execute('select count(*) from %s where id = %%s and expiration_time > now()' % MySQLSession.table_name, (self.id,))
        (count,) = self.cursor.fetchone()
        if count:
            self.cursor.execute('update %s set data = %%s, '
                                'expiration_time = %%s where id = %%s' % MySQLSession.table_name,
                                (pickled_data, expiration_time, self.id))
        else:
            self.cursor.execute('insert into %s (data, expiration_time, id) values (%%s, %%s, %%s)' % MySQLSession.table_name,
                                (pickled_data, expiration_time, self.id))
        self.db.commit()

    def acquire_lock(self):
        logger.debug('acquire_lock %r' % self)
        self.locked = True
        self.cursor.execute('select id from %s where id = %%s for update' % MySQLSession.table_name,
                            (self.id,))
        self.db.commit()

    def release_lock(self):
        logger.debug('release_lock %r' % self)
        self.locked = False
        self.db.commit()

    def clean_up(self):
        logger.debug('clean_up %r' % self)
        self.cursor.execute('delete from %s where expiration_time < now()' % MySQLSession.table_name)
        self.db.commit()

    def _delete(self):
        logger.debug('_delete %r' % self)
        self.cursor.execute('delete from %s where id=%%s' % MySQLSession.table_name, (self.id,))
        self.db.commit()

    def _exists(self):
        # Select session data from table
        self.cursor.execute('select count(*) from %s '
                            'where id = %%s and expiration_time > now()' % MySQLSession.table_name, (self.id,))
        (count,) = self.cursor.fetchone()
        logger.debug('_exists %r (%r)' % (self, bool(count)))
        return bool(count)

    def __del__(self):
        logger.debug('__del__ %r' % self)
        self.db.commit()
        self.db.close()
        self.db = None

    def __repr__(self):
        return '<MySQLSession %r>' % (self.id,)

cherrypy.lib.sessions.MysqlSession = MySQLSession
