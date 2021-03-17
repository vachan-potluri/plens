/**
 * @file BC.h
 * @brief Base class for BC
 */

#include "BC.h"

using namespace BCs;

/**
 * @brief Constructor. Takes and sets the relevant references.
 */
BC::BC(
    const DoFHandler<dim>& dh,
    const std::array<LA::MPI::Vector, 5>& gcv,
    const std::array<LA::MPI::Vector, 9>& gav
):
dof_handler(dh),
g_cvars(gcv),
g_avars(gav),
degree(dh.get_fe().degree),
fdi(dh.get_fe().degree)
{
    form_cell_map();
}



/**
 * @brief Default destructor.
 *
 * This destructor is not made pure virtual to allow testing of the base class BC.
 */
BC::~BC() = default;



/**
 * @brief Populates BC::cell_map_
 *
 * The function loops over active cell iterators and adds the cell index to the map if it is owned.
 */
void BC::form_cell_map()
{
    for(auto cell: dof_handler.active_cell_iterators()){
        if(cell->is_locally_owned()) cell_map_[cell->index()] = cell;
    }
}



/**
 * @brief Returns global dof index based on LocalDoFData object @p ldd
 *
 * Algorithm:
 * - Use cell_map_ to get the cell iterator for cell index `ldd.cell_id`
 * - Get (global) dof indices of this cell
 * - Use FaceDoFInfo @p fdi to get cell-local dof id from `ldd.face_id` and `ldd.face_dof_id`
 * - Use cell dof indices to get global dof from cell-local dof
 */
psize BC::get_global_dof_id(const LocalDoFData &ldd) const
{
    // operator[] of std::map doesn't have a const version, hence use at()
    const DoFHandler<dim>::active_cell_iterator cell = cell_map_.at(ldd.cell_id);
    
    // global dof ids held by this cell
    std::vector<unsigned int> dof_ids(dof_handler.get_fe().dofs_per_cell);
    cell->get_dof_indices(dof_ids);
    
    // get cell local dof id from face local dof id
    const usi cell_local_dof_id = fdi.maps[ldd.face_id][ldd.face_dof_id];
    // get global dof id
    const psize gdof_id = dof_ids[cell_local_dof_id];
    
    return gdof_id;
}



/**
 * @brief Gives state @p s based on @p ldd
 *
 * Internally uses BC::get_global_dof_id()
 */
void BC::get_state(const LocalDoFData &ldd, State &s) const
{
    const psize gdof_id = get_global_dof_id(ldd);
    for(cvar var: cvar_list) s[var] = g_cvars[var][gdof_id];
}



/**
 * @brief Gives avars @p a based on @p ldd
 *
 * Internally uses BC::get_global_dof_id()
 */
void BC::get_avars(const LocalDoFData &ldd, Avars &a) const
{
    const psize gdof_id = get_global_dof_id(ldd);
    for(avar var: avar_list) a[var] = g_avars[var][gdof_id];
}



/**
 * @brief Gives cavars @p ca based on @p ldd
 *
 * Internally uses BC::get_global_dof_id()
 */
void BC::get_cavars(const LocalDoFData &ldd, CAvars &ca) const
{
    const psize gdof_id = get_global_dof_id(ldd);
    
    State& s = ca.get_state();
    Avars& a = ca.get_avars();
    for(cvar var: cvar_list) s[var] = g_cvars[var][gdof_id];
    for(avar var: avar_list) a[var] = g_avars[var][gdof_id];
}



#ifdef DEBUG
void BC::test()
{
    utilities::Testing t("BC", "class");
    utilities::BCTestData bctd(2,2); // refinement and degree
    
    BC bc(bctd.dof_handler, bctd.g_cvars, bctd.g_avars);
    
    {
        t.new_block("testing cavars getter");
        State cons;
        Avars av;
        CAvars ca(&cons,&av);
        DoFHandler<dim>::active_cell_iterator cell =
            bctd.dof_handler.active_cell_iterators().begin();
            // assuming "begin" is owned by this process
        ++cell; // assuming the next iterator is also owned by this process
        LocalDoFData ldd(cell->index(), 1, 3);
        
        psize gdof_id = bc.get_global_dof_id(ldd);
        bc.get_cavars(ldd, ca);
        
        std::cout << "Global dof: " << gdof_id;
        utilities::print_state(ca.get_state());
        utilities::print_avars(ca.get_avars());
    }
}
#endif

