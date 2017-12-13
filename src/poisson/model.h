#pragma once

#include <array>

#include <util/common/plot/plot.h>

#include "plot.h"

namespace model
{

    using points_t = std::vector < plot::point < double > > ;

    struct parameters
    {
        // other params
        double dt, dy, dx;
    };

    inline static parameters make_default_parameters()
    {
        return
        {
        };
    }

    struct plot_data
    {
        util::ptr_t < std::vector < points_t > > data;
        plot::multilist_drawable < points_t > :: ptr_t plot;
        plot::world_t::ptr_t world;
        plot::world_mapper_t world_mapper;
    };

    inline static plot_data make_plot_data()
    {
        plot_data pd;
        pd.data = util::create < std::vector < points_t > > ();
        pd.plot = plot::multilist_drawable < points_t > :: create
        (
            plot::make_data_source(pd.data),
            nullptr, // no point painter
            plot::palette::pen(0xffffff, 2)
        );
        pd.world = plot::world_t::create();
        pd.world_mapper = plot::make_world_mapper(pd.world);
        return pd;
    }

    inline static void adjust(parameters & params,
                              plot_data  & data)
    {
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

    struct chasing_data;

    enum class chasing_dir
    {
        i, j
    };

    struct chasing_coefs
    {
        double A, B, C, D;
    };

    using chasing_fn = std::function < chasing_coefs (size_t i, size_t j, chasing_dir dir,
                                                      const chasing_data & d,
                                                      const parameters & p,
                                                      std::vector < std::vector < double > > & T) > ;

    using material_t = size_t;

    namespace material
    {
        const material_t border_i = 0x1 << 0;
        const material_t border_j = 0x1 << 1;
        const material_t ext      = 0x1 << 4;

        const material_t border   = border_i | border_j;
        const material_t flags    = border;
        const material_t no_flags = ~flags;
    };

    struct chasing_data
    {
        std::vector < std::vector < material_t > > area_map;
        std::vector < std::vector < double > > fn;
        chasing_fn coefs;
        size_t n, m;
    };

    inline material_t get_material_at(const chasing_data & d,
                                      const plot::point < int > & p,
                                      bool preserve_flags = true)
    {
        return ((p.x < 0) || (p.x >= d.n) || (p.y < 0) || (p.y >= d.m))
            ? material::ext
            : (preserve_flags ? d.area_map[p.x][p.y]
                              : (d.area_map[p.x][p.y] & material::no_flags));
    }

    inline std::array < material_t, 2 > get_nearest_materials(const chasing_data & d,
                                                              const plot::point < int > & p,
                                                              chasing_dir dir)
    {
        material_t m1, m2;

        if (dir == chasing_dir::i)
        {
            return {{ get_material_at(d, { p.x - 1, p.y }), get_material_at(d, { p.x + 1, p.y }) }};
        }
        else
        {
            return {{ get_material_at(d, { p.x, p.y - 1 }), get_material_at(d, { p.x, p.y + 1 }) }};
        }
    }

    inline bool is_in_rect(const plot::point < size_t > & p,
                           const plot::rect < size_t > & r)
    {
        return (r.xmin <= p.x) && (p.x <= r.xmax)
            && (r.ymin <= p.y) && (p.y <= r.ymax);
    }

    inline void make_chasing_data(chasing_data & d, const parameters & p)
    {
    }

    // post-apply boundary conditions again
    inline void chasing_bc
    (
        const chasing_data & d,
        const parameters & p,
        std::vector < std::vector < double > > & T
    )
    {
        chasing_coefs c;

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                if (d.area_map[i][j] & material::border_i)
                {
                    if ((get_material_at(d, { (int) i - 1, (int) j })
                        | get_material_at(d, { (int) i + 1, (int) j })) & material::ext)
                    {
                        c = d.coefs(i, j, chasing_dir::i, d, p, T);
                        T[i][j] = c.D;
                        if (!(get_material_at(d, { (int) i - 1, (int) j }) & material::ext))
                            T[i][j] -= c.A * T[i - 1][j];
                        if (!(get_material_at(d, { (int) i + 1, (int) j }) & material::ext))
                            T[i][j] -= c.B * T[i + 1][j];
                        T[i][j] /= c.C;
                    }
                    else if ((get_material_at(d, { (int) i, (int) j - 1 })
                        | get_material_at(d, { (int) i, (int) j + 1 })) & material::ext)
                    {
                        c = d.coefs(i, j, chasing_dir::j, d, p, T);
                        T[i][j] = c.D;
                        if (!(get_material_at(d, { (int) i, (int) j - 1 }) & material::ext))
                            T[i][j] -= c.A * T[i][j - 1];
                        if (!(get_material_at(d, { (int) i, (int) j + 1 }) & material::ext))
                            T[i][j] -= c.B * T[i][j + 1];
                        T[i][j] /= c.C;
                    }
                }
            }
        }
    }

    inline void chasing_solve
    (
        const chasing_data & d,
        const parameters & p,
        std::vector < std::vector < double > > & T
    )
    {
        std::vector < std::vector < double > > a(d.n + 1), b(d.n + 1);

        for (size_t i = 0; i < d.n + 1; ++i)
        {
            a[i].resize(d.m + 1); b[i].resize(d.m + 1);
        }

        chasing_coefs c;

        // i-chasing

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                c = d.coefs(i, j, chasing_dir::i, d, p, T);
                a[i + 1][j] = - c.B / (c.C + c.A * a[i][j]);
                b[i + 1][j] = (c.D - c.A * b[i][j]) / (c.C + c.A * a[i][j]);
            }
        }

        for (size_t j = 0; j < d.m; ++j)
        {
            size_t s = d.n - 1, e = 0;
            while (get_material_at(d, { (int) s, (int) j }) & material::ext) --s;
            while (get_material_at(d, { (int) e, (int) j }) & material::ext) ++e;
            for (size_t i = s + 2; i-- > e + 1;)
            {
                if (get_material_at(d, { (int) i, (int) j }) & material::ext)
                {
                    T[i - 1][j] = b[i][j];
                }
                else
                {
                    T[i - 1][j] = a[i][j] * T[i][j] + b[i][j];
                }
            }
        }

        // boundary condition chasing

        chasing_bc(d, p, T);

        // j-chasing

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                c = d.coefs(i, j, chasing_dir::j, d, p, T);
                a[i][j + 1] = - c.B / (c.C + c.A * a[i][j]);
                b[i][j + 1] = (c.D - c.A * b[i][j]) / (c.C + c.A * a[i][j]);
            }
        }

        for (size_t i = 0; i < d.n; ++i)
        {
            size_t s = d.m - 1, e = 0;
            while (get_material_at(d, { (int) i, (int) s }) & material::ext) --s;
            while (get_material_at(d, { (int) i, (int) e }) & material::ext) ++e;
            for (size_t j = s + 2; j-- > e + 1;)
            {
                if (get_material_at(d, { (int) i, (int) j }) & material::ext)
                {
                    T[i][j - 1] = b[i][j];
                }
                else
                {
                    T[i][j - 1] = a[i][j] * T[i][j] + b[i][j];
                }
            }
        }

        // boundary condition chasing

        chasing_bc(d, p, T);
    }

    using stencil_fn = std::function < bool (int i, int j) > ;

    inline static stencil_fn make_simple_stencil(size_t n, size_t m)
    {
        return [=] (int i, int j) { return (i >= 0) && (j >= 0) && (i < n) && (j < m); };
    }

    inline static stencil_fn make_material_based_stencil(const chasing_data & d)
    {
        return [&] (int i, int j) { return !(get_material_at(d, { i, j }) & material::ext); };
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

    inline static plot::painter_t make_system_painter(const parameters & params,
                                                      const chasing_data & d,
                                                      const std::vector < std::vector < double > > & T,
                                                      double & max_heat_value,
                                                      bool & draw_heat_map)
    {
        return [&] (CDC & dc, const plot::viewport & vp)
        {
            if (d.area_map.empty()) return;
        };
    }
}
