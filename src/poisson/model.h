#pragma once

#include <array>
#include <map>

#include <util/common/plot/plot.h>

#include "plot.h"

namespace model
{

    using points_t = std::vector < plot::point < double > > ;

    struct parameters
    {
        // system params
        double a, b, k;

        // capacitor params
        double d, angle;

        // other params
        double dt, dy, dx;

        // material params
        double eps, q1, q2;
    };

    inline static parameters make_default_parameters()
    {
        return
        {
            // system params
            50, 50, 2,

            // capacitor params
            5, 30,

            // other params
            0.3, 1, 1,

            // material params
            1, 10, -10
        };
    }

    struct plot_data
    {
        util::ptr_t < std::vector < points_t > > data;
        plot::multilist_drawable < points_t > :: ptr_t plot;
        plot::world_t::ptr_t world;
        plot::world_mapper_t world_mapper;
    };

    inline static plot_data make_plot_data(COLORREF color = 0xffffff)
    {
        plot_data pd;
        pd.data = util::create < std::vector < points_t > > ();
        pd.plot = plot::multilist_drawable < points_t > :: create
        (
            plot::make_data_source(pd.data),
            nullptr, // no point painter
            plot::palette::pen(color, 2)
        );
        pd.world = plot::world_t::create();
        pd.world_mapper = plot::make_world_mapper(pd.world);
        return pd;
    }

    inline static void adjust(parameters & params,
                              plot_data  & data)
    {
        data.world->xmin = - (data.world->xmax = params.a * params.k);
        data.world->ymin = - (data.world->ymax = params.b * params.k);
    }

    inline static plot::drawable::ptr_t make_root_drawable
    (
        const plot_data & data,
        std::vector < plot::drawable::ptr_t > layers
    )
    {
        using namespace plot;

        layers.push_back(data.plot);

        return viewporter::create(
            tick_drawable::create(
                layer_drawable::create(layers),
                const_n_tick_factory<axe::x>::create(
                    make_simple_tick_formatter(2, 5),
                    0,
                    5
                ),
                const_n_tick_factory<axe::y>::create(
                    make_simple_tick_formatter(2, 5),
                    0,
                    5
                ),
                palette::pen(RGB(80, 80, 80)),
                RGB(200, 200, 200)
            ),
            make_viewport_mapper(data.world_mapper)
        );
    }

    inline static RECT xyxy(const plot::viewport & vp, const plot::rect < double > & r)
    {
        return
        {
            vp.world_to_screen().x(r.xmin),
            vp.world_to_screen().y(r.ymin),
            vp.world_to_screen().x(r.xmax),
            vp.world_to_screen().y(r.ymax)
        };
    }

    struct relax_data;

    using material_t = size_t;

    namespace material
    {
        const material_t border_i = 0x1 << 0;
        const material_t border_j = 0x1 << 1;
        const material_t ext      = 0x1 << 4;
        const material_t metal    = 0x1 << 5;
        const material_t dielectr = 0x1 << 6;
        const material_t cap1     = 0x1 << 7;
        const material_t cap2     = 0x1 << 8;

        const material_t border   = border_i | border_j;
        const material_t flags    = border;
        const material_t no_flags = ~flags;
    };

    struct relax_data
    {
        std::vector < std::vector < material_t > > area_map;
        std::vector < std::vector < double > > fn;
        size_t n, m;
    };

    inline material_t get_material_at(const relax_data & d,
                                      const plot::point < int > & p,
                                      bool preserve_flags = true)
    {
        return ((p.x < 0) || (p.x >= d.n) || (p.y < 0) || (p.y >= d.m))
            ? material::ext
            : (preserve_flags ? d.area_map[p.x][p.y]
                              : (d.area_map[p.x][p.y] & material::no_flags));
    }

    inline bool is_in_rect(const plot::point < size_t > & p,
                           const plot::rect < size_t > & r)
    {
        return (r.xmin <= p.x) && (p.x <= r.xmax)
            && (r.ymin <= p.y) && (p.y <= r.ymax);
    }

    inline void make_relax_data(relax_data & d, const parameters & p)
    {
        d.n = (size_t) std::ceil((p.b * p.k) / p.dy) * 2;
        d.m = (size_t) std::ceil((p.a * p.k) / p.dx) * 2;

        d.area_map.clear();
        d.fn.clear();

        d.area_map.resize(d.n, std::vector < material_t > (d.m));
        d.fn.resize(d.n, std::vector < double > (d.m));

        size_t Y_n  = (size_t) std::ceil((p.b * p.k) / p.dy);
        size_t X_m  = (size_t) std::ceil((p.a * p.k) / p.dx);
        size_t b_n  = (size_t) std::ceil(p.b / p.dy);
        size_t a_m  = (size_t) std::ceil(p.a / p.dx);
        size_t d_n  = (size_t) std::ceil(p.d / p.dy);
        size_t d_m  = (size_t) std::ceil(p.d / p.dx);

        // set up materials and sketch out borders

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                bool is_border           = (i == 0) || ((i + 1) == Y_n * 2) || (j == 0) || ((j + 1) == X_m * 2);
                bool is_capacitor        = is_in_rect({ i, j }, { Y_n - b_n, Y_n + b_n, X_m - a_m, X_m - a_m + d_m })
                                        || is_in_rect({ i, j }, { Y_n - b_n, Y_n + b_n, X_m + a_m - d_m, X_m + a_m });
                bool is_capacitor_border = is_capacitor && (
                                                !is_in_rect({ i, j }, { Y_n - b_n + 1, Y_n + b_n - 1, X_m - a_m + 1, X_m - a_m + d_m - 1 })
                                             && !is_in_rect({ i, j }, { Y_n - b_n + 1, Y_n + b_n - 1, X_m + a_m - d_m + 1, X_m + a_m - 1 })
                                         );

                if (is_border || is_capacitor_border)
                    d.area_map[i][j]  = material::border;
                if (is_capacitor || is_border)
                    d.area_map[i][j] |= material::metal;
                else
                    d.area_map[i][j]  = material::dielectr;

                if (is_capacitor && (j < Y_n))
                {
                    d.fn[i][j] = p.q1;
                    d.area_map[i][j] |= material::cap1;
                }
                if (is_capacitor && (j > Y_n))
                {
                    d.fn[i][j] = p.q2;
                    d.area_map[i][j] |= material::cap2;
                }
            }
        }

        // detect border orientations

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                if (get_material_at(d, { (int) i, (int) j }) & material::border)
                {
                    d.area_map[i][j] &= ~material::border;
                    if ((get_material_at(d, { (int) i - 1, (int) j })
                        | get_material_at(d, { (int) i + 1, (int) j })) & material::border)
                    {
                        d.area_map[i][j] |= material::border_j;
                    }
                    if ((get_material_at(d, { (int) i, (int) j - 1 })
                        | get_material_at(d, { (int) i, (int) j + 1 })) & material::border)
                    {
                        d.area_map[i][j] |= material::border_i;
                    }
                }
            }
        }
    }

    inline void relax_solve
    (
        const relax_data & d,
        const parameters & p,
        std::vector < std::vector < double > > & T
    )
    {
        for (size_t l = 0; l < 100; ++l)
        {
            for (size_t i = 0; i < d.n; ++i)
            {
                for (size_t j = 0; j < d.m; ++j)
                {
                    if (d.area_map[i][j] & material::border)
                    {
                             if (d.area_map[i][j] & material::cap1) T[i][j] = p.q1;
                        else if (d.area_map[i][j] & material::cap2) T[i][j] = p.q2;
                        else T[i][j] = 0;
                    }
                    else if (d.area_map[i][j] & material::metal) continue;
                    else
                    {
                        double Tij = 0;
                        if (!(get_material_at(d, { (int) i - 1, (int) j }) & material::ext))
                            Tij += T[i - 1][j] / p.dy / p.dy;
                        if (!(get_material_at(d, { (int) i + 1, (int) j }) & material::ext))
                            Tij += T[i + 1][j] / p.dy / p.dy;
                        if (!(get_material_at(d, { (int) i, (int) j - 1 }) & material::ext))
                            Tij += T[i][j - 1] / p.dx / p.dx;
                        if (!(get_material_at(d, { (int) i, (int) j + 1 }) & material::ext))
                            Tij += T[i][j + 1] / p.dx / p.dx;
                        Tij -= d.fn[i][j] / p.eps;
                        Tij /= 2 * (1. / p.dx / p.dx + 1. / p.dy / p.dy);
                        T[i][j] = Tij;
                    }
                }
            }
        }
    }

    using stencil_fn = std::function < bool (int i, int j) > ;

    inline static stencil_fn make_simple_stencil(size_t n, size_t m)
    {
        return [=] (int i, int j) { return (i >= 0) && (j >= 0) && (i < n) && (j < m); };
    }

    inline static stencil_fn make_material_based_stencil(const relax_data & d)
    {
        return [&] (int i, int j) { return !(get_material_at(d, { i, j }) & (material::ext | material::metal)); };
    }

    void find_isolines
    (
        const std::vector < std::vector < double > > & T,
        double dT,
        std::vector < std::vector < plot::point < double > > > & out,
        size_t n, size_t m,
        const parameters & p,
        stencil_fn stencil,
        size_t max_isolines = 100,
        size_t max_points_in_stack = 100000
    );

    void find_field_lines
    (
        const std::vector < std::vector < double > > & T,
        std::vector < std::vector < plot::point < double > > > & out,
        size_t n, size_t m,
        const parameters & p,
        stencil_fn stencil,
        std::vector < plot::point < double > > hint
    );

    inline static plot::painter_t make_system_painter(const parameters & params,
                                                      const relax_data & d,
                                                      const std::vector < std::vector < double > > & T)
    {
        return [&] (CDC & dc, const plot::viewport & vp)
        {
            if (d.area_map.empty()) return;

            auto metal_brush  = plot::palette::brush(RGB(100, 100, 100));
            auto border_brush = plot::palette::brush(RGB(155, 155, 155));

            RECT r;

            material_t m;

            for (size_t i = 0; i < d.n; ++i)
            {
                for (size_t j = 0; j < d.m; ++j)
                {
                    m = get_material_at(d, { (int) i, (int) j });

                    if (m & material::ext) continue;

                    if (m & (material::border | material::metal))
                    {
                        r = xyxy(vp,
                        {
                            ((double)j - 0.5) * params.dx - params.a * params.k,
                            ((double)j + 0.5) * params.dx - params.a * params.k,
                            ((double)i - 0.5) * params.dy - params.b * params.k,
                            ((double)i + 0.5) * params.dy - params.b * params.k
                        });

                        if (m & material::border)
                        {
                            dc.FillRect(&r, border_brush.get());
                        }
                        else if (m & material::metal)
                        {
                            dc.FillRect(&r, metal_brush.get());
                        }
                    }
                }
            }
        };
    }
}
