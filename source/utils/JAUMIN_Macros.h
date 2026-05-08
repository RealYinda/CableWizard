#ifndef included_JAUMIN_Macros
#define included_JAUMIN_Macros

// NOTE:
//      The following macros are used to simplify coding applications
//      on JAUMIN framework. However, they are not part of the framework.
//      If you decide to use these macros, follow the examples, just guess their
//      meanings from the examples, don't try to read their definition.
//      If you are not familiar with JAUMIN framework, DON'T use these macros.

#define DECLARE_VARIABLE(var, entity, type, group, depth)  \
  tbox::Pointer<pdat::entity##Variable<NDIM, type> > var = \
      new pdat::entity##Variable<NDIM, type>(#var, group, depth);

#define DECLARE_MATVEC_VARIABLE(var, matvec, type, dofinfo) \
  tbox::Pointer<pdat::matvec##Variable<NDIM, type> > var =  \
      new pdat::matvec##Variable<NDIM, type>(#var, dofinfo);

#define REGISTER_VARIABLE(id, var, context, ghost)                     \
  {                                                                    \
    hier::VariableDatabase<NDIM>* variable_db =                        \
        hier::VariableDatabase<NDIM>::getDatabase();                   \
    tbox::Pointer<hier::VariableContext> context =                     \
        variable_db->getContext(#context);                             \
    id = variable_db->registerVariableAndContext(var, context, ghost); \
  }

#define DECLARE_ADJACENCY(patch, entity1, entity2, Entity1, Entity2)     \
  tbox::Array<int> entity1##_##entity2##_ext, entity1##_##entity2##_idx; \
  patch.getPatchTopology()->get##Entity1##Adjacency##Entity2##s(         \
      entity1##_##entity2##_ext, entity1##_##entity2##_idx);

#define GET_PATCH_DATA(patch, data, id, entity, type)          \
  tbox::Pointer<pdat::entity##Data<NDIM, type> > data =        \
      patch.getPatchData(id);                                  \
  if (data.isNull()) {                                         \
    TBOX_ERROR("Mismatch data type for \"" << #id << "\".\n"); \
  }

#define GET_COORD_DATA(patch, data, entity)               \
  tbox::Pointer<pdat::entity##Data<NDIM, double> > data = \
      patch.getPatchGeometry()->get##entity##Coordinates();

#define HAS_ENTITY_SET(patch, set_id, ETYPE, ghost) \
  patch.getPatchGeometry()->hasEntitySet(           \
      set_id, hier::EntityUtilities::ETYPE,         \
      patch.getNumberOfEntities(hier::EntityUtilities::ETYPE, ghost))

#define DECLARE_ENTITY_SET(patch, set, set_id, ETYPE, ghost) \
  const tbox::Array<int>& set =                              \
      patch.getPatchGeometry()->getEntityIndicesInSet(       \
          set_id, hier::EntityUtilities::ETYPE,              \
          patch.getNumberOfEntities(hier::EntityUtilities::ETYPE, ghost));
          
#define DOT_PRODUCT(v1, v2) \
  v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];

#define TIMES(v1, v2, v3) \
    v3[0] = v1*v2[0]; \
    v3[1] = v1*v2[1]; \
    v3[2] = v1*v2[2];

#define CROSS_PRODUCT(a, b, n)  \
  n[0] = a[1] * b[2] - a[2] * b[1]; \
  n[1] = a[2] * b[0] - a[0] * b[2]; \
  n[2] = a[0] * b[1] - a[1] * b[0]; 
  
#define NORMALIZATION(n)  \
  { \
  double dn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];  \
  dn = 1.0 / (dn == 0.0 ? 1.0 : sqrt(dn)); \
  n[0] *= dn; \
  n[1] *= dn; \
  n[2] *= dn; \
  }
  
#define POUT_3D_VECTOR(v1)  \
  tbox::pout << #v1 << " : " << v1[0] << ", " << v1[1] << ", " <<v1[2] <<endl;

#define USERDEFINEGRADIENT(gradient, p0, p1, vol)          \
  {                                              \
    (gradient)[0] = ((p0)[1] - (p1)[1]) / (vol); \
    (gradient)[1] = ((p1)[0] - (p0)[0]) / (vol); \
  }
  
#endif  // included_JAUMIN_Macros
