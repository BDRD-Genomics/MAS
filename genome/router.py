"""
This is the router that will send the genome models to the genome database
"""
# class GenomesRouter(object):
#     #Suggest the database that should be used for read operations
#     def db_for_read(self,model,**hints):
#         if model._meta.app_label == 'genome':
#             #genome db
#             return 'genome'
#         return None
#
#     #Suggest the database that should be used for writes
#     def db_for_write(self,model,**hints):
#         if model._meta.app_label == 'genome':
#             return 'genome'
#         return None
#
#     #Return True if a relation between obj1 and obj2 should be allowed, False if the relation should be prevented,
#     # or None if the router has no opinion.
#     def allow_relation(self,obj1,obj2,**hints):
#         # if obj1._meta.app_label == 'genome' or \
#         #         obj2._meta.app_label == 'genome':
#         #     return True
#         if obj1._meta.app_label == 'genome' and obj2._meta.app_label == 'genome':
#             return True
#         elif 'genome' not in [obj1._meta.app_label,obj2._meta.app_label]:
#             return None
#
#         return False
#
#     #Determine if the migration operation is allowed to run on the database with alias db.
#     # Return True if the operation should run, False if it shouldn’t run, or None if the router has no opinion.
#     def allow_migrate(self,db,app_label,model_name=None,**hints):
#         if app_label == 'genome':
#             return db == 'genome'
#         elif db == 'genome':
#             return False
#         return None

# """
# This is the router that will send the amd_models to the amd_genome database
# """
# class AmdRouter(object):
#     #Suggest the database that should be used for read operations
#     def db_for_read(self,model,**hints):
#         if model._meta.app_label == 'amd_genome':
#             #genome db
#             return 'amd_genome'
#         return None
#
#     #Suggest the database that should be used for writes
#     def db_for_write(self,model,**hints):
#         if model._meta.app_label == 'amd_genome':
#             return 'amd_genome'
#         return None
#
#     #Return True if a relation between obj1 and obj2 should be allowed, False if the relation should be prevented,
#     # or None if the router has no opinion.
#     def allow_relation(self,obj1,obj2,**hints):
#         # if obj1._meta.app_label == 'genome' or \
#         #         obj2._meta.app_label == 'genome':
#         #     return True
#         if obj1._meta.app_label == 'amd_genome' and obj2._meta.app_label == 'amd_genome':
#             return True
#         elif 'amd_genome' not in [obj1._meta.app_label,obj2._meta.app_label]:
#             return None
#
#         return False
#
#     #Determine if the migration operation is allowed to run on the database with alias db.
#     # Return True if the operation should run, False if it shouldn’t run, or None if the router has no opinion.
#     def allow_migrate(self,db,app_label,model_name=None,**hints):
#         if app_label == 'amd_genome':
#             return db == 'amd_genome'
#         elif db == 'amd_genome':
#             return False
#         return None