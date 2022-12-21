static const char help[] = "Simple libCEED test to calculate finite volume residual";

#include <petscdmceed.h>
#include <petscdmplexceed.h>
#include <petscfeceed.h>
#include <petscfvceed.h>
#include <petscdmplex.h>
#include <petscds.h>

typedef struct {
  CeedQFunctionUser setupgeo, apply;
  const char       *setupgeofname, *applyfname;
} AppCtx;

typedef struct {
  CeedQFunction qf_apply;
  CeedOperator  op_apply;
  CeedVector    qdata, uceed, vceed;
} CeedData;

CEED_QFUNCTION(SetupGeo)(void *ctx, const CeedInt Q, const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *J = in[1], *w = in[2];
  CeedScalar       *qdata = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; ++i)
  {
    // Read dxxdX Jacobian entries, stored as [[0 3], [1 4], [2 5]]
    const CeedScalar dxxdX[3][2] = {
      {J[i + Q * 0], J[i + Q * 3]},
      {J[i + Q * 1], J[i + Q * 4]},
      {J[i + Q * 2], J[i + Q * 5]}
    };
    // Modulus of dxxdX column vectors
    const CeedScalar modg1 = PetscSqrtReal(dxxdX[0][0] * dxxdX[0][0] + dxxdX[1][0] * dxxdX[1][0] + dxxdX[2][0] * dxxdX[2][0]);
    const CeedScalar modg2 = PetscSqrtReal(dxxdX[0][1] * dxxdX[0][1] + dxxdX[1][1] * dxxdX[1][1] + dxxdX[2][1] * dxxdX[2][1]);
    // Use normalized column vectors of dxxdX as rows of dxdxx
    const CeedScalar dxdxx[2][3] = {
      {dxxdX[0][0] / modg1, dxxdX[1][0] / modg1, dxxdX[2][0] / modg1},
      {dxxdX[0][1] / modg2, dxxdX[1][1] / modg2, dxxdX[2][1] / modg2}
    };

    CeedScalar dxdX[2][2];
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k) {
        dxdX[j][k] = 0;
        for (int l = 0; l < 3; ++l) dxdX[j][k] += dxdxx[j][l] * dxxdX[l][k];
      }
    qdata[i + Q * 0] = (dxdX[0][0] * dxdX[1][1] - dxdX[1][0] * dxdX[0][1]) * w[i]; /* det J * weight */
  }
  return 0;
}

CEED_QFUNCTION(Advection)(void *ctx, const CeedInt Q, const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *u = in[0], *qdata = in[1];
  CeedScalar       *v = out[0];

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; ++i) v[i] = qdata[i] * u[i];

  return 0;
}

static PetscErrorCode CeedDataDestroy(CeedData *data)
{
  PetscFunctionBeginUser;
  PetscCall(CeedVectorDestroy(&data->qdata));
  PetscCall(CeedVectorDestroy(&data->uceed));
  PetscCall(CeedVectorDestroy(&data->vceed));
  PetscCall(CeedQFunctionDestroy(&data->qf_apply));
  PetscCall(CeedOperatorDestroy(&data->op_apply));
  PetscFunctionReturn(0);
}

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *ctx)
{
  PetscFunctionBeginUser;
  PetscOptionsBegin(comm, "", "libCEED Test Options", "DMPLEX");
  PetscOptionsEnd();
  ctx->setupgeo      = SetupGeo;
  ctx->setupgeofname = SetupGeo_loc;
  ctx->apply         = Advection;
  ctx->applyfname    = Advection_loc;
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateMesh(MPI_Comm comm, AppCtx *ctx, DM *dm)
{
  PetscFunctionBegin;
  PetscCall(DMCreate(comm, dm));
  PetscCall(DMSetType(*dm, DMPLEX));
  PetscCall(DMSetFromOptions(*dm));
  PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));
#ifdef PETSC_HAVE_LIBCEED
  {
    Ceed        ceed;
    const char *usedresource;

    PetscCall(DMGetCeed(*dm, &ceed));
    PetscCall(CeedGetResource(ceed, &usedresource));
    PetscCall(PetscPrintf(PetscObjectComm((PetscObject)*dm), "libCEED Backend: %s\n", usedresource));
  }
#endif
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupDiscretization(DM dm)
{
  PetscFV  fv;
  PetscInt dim;

  PetscFunctionBeginUser;
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(PetscFVCreate(PetscObjectComm((PetscObject)dm), &fv));
  PetscCall(PetscFVSetFromOptions(fv));
  PetscCall(PetscFVSetNumComponents(fv, 1));
  PetscCall(PetscFVSetSpatialDimension(fv, dim));
  PetscCall(PetscObjectSetName((PetscObject)fv, "concentration"));
  PetscCall(DMAddField(dm, NULL, (PetscObject)fv));
  PetscCall(PetscFVDestroy(&fv));
  PetscCall(DMCreateDS(dm));
  PetscFunctionReturn(0);
}

static PetscErrorCode LibCeedSetupByDegree(DM dm, AppCtx *ctx, CeedData *data)
{
  PetscDS             ds;
  PetscFV             fv;
  PetscFE             cfe;
  Ceed                ceed;
  CeedElemRestriction Erestrictx, Erestrictneg, Erestrictpos, Erestrictq;
  CeedQFunction       qf_setupgeo;
  CeedOperator        op_setupgeo;
  CeedVector          xcoord;
  CeedBasis           basisu, basisx;
  CeedInt             Nqdata = 1;
  CeedInt             nqpts, nqptsx;
  DM                  cdm;
  Vec                 coords;
  const PetscScalar  *coordArray;
  PetscInt            dim, cdim, cStart, cEnd, Ncell;

  PetscFunctionBeginUser;
  PetscCall(DMGetCeed(dm, &ceed));
  PetscCall(DMGetDimension(dm, &dim));
  PetscCall(DMGetCoordinateDim(dm, &cdim));
  PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
  Ncell = cEnd - cStart;
  // CEED bases
  PetscCall(DMGetDS(dm, &ds));
  PetscCall(PetscDSGetDiscretization(ds, 0, (PetscObject *)&fv));
  PetscCall(PetscFVGetCeedBasis(fv, &basisu));
  PetscCall(DMGetCoordinateDM(dm, &cdm));
  PetscCall(DMGetDS(cdm, &ds));
  PetscCall(PetscDSGetDiscretization(ds, 0, (PetscObject *)&cfe));
  PetscCall(PetscFEGetCeedBasis(cfe, &basisx));

  PetscCall(DMPlexGetCeedRestriction(cdm, NULL, 0, 0, 0, &Erestrictx));
  PetscCall(DMPlexGetCeedFaceRestriction(dm, NULL, 0, 0, 0, &Erestrictneg, &Erestrictpos));
  PetscCall(CeedBasisGetNumQuadraturePoints(basisu, &nqpts));
  PetscCall(CeedBasisGetNumQuadraturePoints(basisx, &nqptsx));
  PetscCheck(nqptsx == nqpts, PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Number of qpoints for u %" PetscInt_FMT " != %" PetscInt_FMT " Number of qpoints for x", nqpts, nqptsx);
  PetscCall(CeedElemRestrictionCreateStrided(ceed, Ncell, nqpts, Nqdata, Nqdata * Ncell * nqpts, CEED_STRIDES_BACKEND, &Erestrictq));

  PetscCall(DMGetCoordinatesLocal(dm, &coords));
  PetscCall(VecGetArrayRead(coords, &coordArray));
  PetscCall(CeedElemRestrictionCreateVector(Erestrictx, &xcoord, NULL));
  PetscCall(CeedVectorSetArray(xcoord, CEED_MEM_HOST, CEED_COPY_VALUES, (PetscScalar *)coordArray));
  PetscCall(VecRestoreArrayRead(coords, &coordArray));

  // Create the vectors that will be needed in setup and apply
  PetscCall(CeedElemRestrictionCreateVector(Erestrictneg, &data->uceed, NULL));
  PetscCall(CeedElemRestrictionCreateVector(Erestrictneg, &data->vceed, NULL));
  PetscCall(CeedElemRestrictionCreateVector(Erestrictq, &data->qdata, NULL));

  // Create the Q-function that builds the operator (i.e. computes its quadrature data) and set its context data
  PetscCall(CeedQFunctionCreateInterior(ceed, 1, ctx->setupgeo, ctx->setupgeofname, &qf_setupgeo));
  PetscCall(CeedQFunctionAddInput(qf_setupgeo, "x", cdim, CEED_EVAL_INTERP));
  PetscCall(CeedQFunctionAddInput(qf_setupgeo, "dx", cdim * dim, CEED_EVAL_GRAD));
  PetscCall(CeedQFunctionAddInput(qf_setupgeo, "weight", 1, CEED_EVAL_WEIGHT));
  PetscCall(CeedQFunctionAddOutput(qf_setupgeo, "qdata", Nqdata, CEED_EVAL_NONE));

  // Set up the residual operator
  PetscCall(CeedQFunctionCreateInterior(ceed, 1, ctx->apply, ctx->applyfname, &data->qf_apply));
  PetscCall(CeedQFunctionAddInput(data->qf_apply, "u", 1, CEED_EVAL_INTERP));
  PetscCall(CeedQFunctionAddInput(data->qf_apply, "qdata", Nqdata, CEED_EVAL_NONE));
  PetscCall(CeedQFunctionAddOutput(data->qf_apply, "v", 1, CEED_EVAL_INTERP));

  // Create the operator that builds the quadrature data for the operator
  PetscCall(CeedOperatorCreate(ceed, qf_setupgeo, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &op_setupgeo));
  PetscCall(CeedOperatorSetField(op_setupgeo, "x", Erestrictx, basisx, CEED_VECTOR_ACTIVE));
  PetscCall(CeedOperatorSetField(op_setupgeo, "dx", Erestrictx, basisx, CEED_VECTOR_ACTIVE));
  PetscCall(CeedOperatorSetField(op_setupgeo, "weight", CEED_ELEMRESTRICTION_NONE, basisx, CEED_VECTOR_NONE));
  PetscCall(CeedOperatorSetField(op_setupgeo, "qdata", Erestrictq, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE));

  // Create the residual operator
  PetscCall(CeedOperatorCreate(ceed, data->qf_apply, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &data->op_apply));
  PetscCall(CeedOperatorSetField(data->op_apply, "uneg", Erestrictneg, basisu, CEED_VECTOR_ACTIVE));
  PetscCall(CeedOperatorSetField(data->op_apply, "upos", Erestrictpos, basisu, CEED_VECTOR_ACTIVE));
  PetscCall(CeedOperatorSetField(data->op_apply, "qdata", Erestrictq, CEED_BASIS_COLLOCATED, data->qdata));
  PetscCall(CeedOperatorSetField(data->op_apply, "fneg", Erestrictneg, basisu, CEED_VECTOR_ACTIVE));
  PetscCall(CeedOperatorSetField(data->op_apply, "fpos", Erestrictpos, basisu, CEED_VECTOR_ACTIVE));

  // Setup qdata
  PetscCall(CeedOperatorApply(op_setupgeo, xcoord, data->qdata, CEED_REQUEST_IMMEDIATE));

  // JED: Who destroys the restrictions? It causes a SEGV to destroy them here
  PetscCall(CeedElemRestrictionDestroy(&Erestrictq));
  PetscCall(CeedQFunctionDestroy(&qf_setupgeo));
  PetscCall(CeedOperatorDestroy(&op_setupgeo));
  PetscCall(CeedVectorDestroy(&xcoord));
  PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
  MPI_Comm     comm;
  DM           dm;
  AppCtx       ctx;
  Vec          U, Uloc, F, Floc;
  PetscScalar *u, *f;
  CeedData     ceeddata;

  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  comm = PETSC_COMM_WORLD;
  PetscCall(ProcessOptions(comm, &ctx));
  PetscCall(CreateMesh(comm, &ctx, &dm));
  PetscCall(SetupDiscretization(dm));

  PetscCall(LibCeedSetupByDegree(dm, &ctx, &ceeddata));

  PetscCall(DMCreateGlobalVector(dm, &U));
  PetscCall(DMCreateLocalVector(dm, &Uloc));
  PetscCall(VecDuplicate(U, &F));
  PetscCall(VecDuplicate(Uloc, &Floc));

  PetscCall(VecSet(Uloc, 1.));
  PetscCall(VecZeroEntries(F));
  PetscCall(VecZeroEntries(Floc));
  PetscCall(VecGetArray(Uloc, &u));
  PetscCall(VecGetArray(Floc, &f));
  PetscCall(CeedVectorSetArray(ceeddata.uceed, CEED_MEM_HOST, CEED_USE_POINTER, u));
  PetscCall(CeedVectorSetArray(ceeddata.vceed, CEED_MEM_HOST, CEED_USE_POINTER, f));
  PetscCall(CeedOperatorApply(ceeddata.op_apply, ceeddata.uceed, ceeddata.vceed, CEED_REQUEST_IMMEDIATE));
  PetscCall(CeedVectorTakeArray(ceeddata.vceed, CEED_MEM_HOST, NULL));
  PetscCall(VecRestoreArray(Floc, &f));
  PetscCall(VecRestoreArray(Uloc, &u));
  PetscCall(DMLocalToGlobalBegin(dm, Floc, ADD_VALUES, F));
  PetscCall(DMLocalToGlobalEnd(dm, Floc, ADD_VALUES, F));

  PetscCall(CeedDataDestroy(&ceeddata));
  PetscCall(VecDestroy(&U));
  PetscCall(VecDestroy(&Uloc));
  PetscCall(VecDestroy(&F));
  PetscCall(VecDestroy(&Floc));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}

/*TEST

  build:
    requires: libceed

  testset:
    args: -dm_plex_simplex 0 -dm_view -dm_petscds_view

    test:
      suffix: square_0
      args: -dm_refine 2

TEST*/
