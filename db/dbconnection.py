from contextlib import contextmanager
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

class DBconnector:
    def __init__(self):
        self.engine = create_engine('mysql+pymysql://connection_credential')
        self.Session = sessionmaker()
        self.Session.configure(bind=self.engine)
        self.session = ''

    @contextmanager
    # create a decorator for db connection, always use it as: with create_session() as xxx
    def create_session(self):
        self.session = self.Session()
        try:
            yield self.session
            self.session.commit()
        except:
            self.session.rollback()
            raise
        finally:
            self.session.close()

   
def main():
    dbc = DBconnector()


if __name__ == "__main__":
    main()
