class DatabaseRouter:
    def db_for_read(self, model, **hints):
        if hasattr(model, 'db'):
            return model.db
        else:
            return 'default'

    def db_for_write(self, model, **hints):
        if hasattr(model, 'db'):
            return model.db
        else:
            return 'default'

    def allow_relation(self, obj1, obj2, **hints):
        if hasattr(obj1, 'db') and hasattr(obj2, 'db') and obj1.db == obj2.db:
            return True

        elif not hasattr(obj1, 'db') and not hasattr(obj2, 'db'):
            return True

        else:
            return False

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        if db == 'AMD':
            return False
        return None
