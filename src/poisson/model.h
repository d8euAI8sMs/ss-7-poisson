#pragma once

#include <array>

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
        double d;

        // Faraday cage params
        double d_c, x_c, y_c, a_c, b_c;

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
            5,

            // Faraday cage params
            5, 0, 0, 20, 20,

            // other params
            0.3, 1, 1,

            // material params
            1, 1, 1
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
        const material_t metal    = 0x1 << 5;
        const material_t dielectr = 0x1 << 6;
        const material_t cap1     = 0x1 << 7;
        const material_t cap2     = 0x1 << 8;

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
        d.n = (size_t) std::ceil((p.b * p.k) / p.dy) * 2;
        d.m = (size_t) std::ceil((p.a * p.k) / p.dx) * 2;

        d.area_map.clear();
        d.fn.clear();

        d.area_map.resize(d.n, std::vector < material_t > (d.m));
        d.fn.resize(d.n, std::vector < double > (d.m));

        size_t Y_n  = (size_t) std::ceil((p.b * p.k) / p.dy);
        size_t X_m  = (size_t) std::ceil((p.a * p.k) / p.dx);
        size_t yc_n = (size_t) std::ceil(p.y_c / p.dy);
        size_t xc_m = (size_t) std::ceil(p.x_c / p.dx);
        size_t b_n  = (size_t) std::ceil(p.b / p.dy);
        size_t a_m  = (size_t) std::ceil(p.a / p.dx);
        size_t bc_n  = (size_t) std::ceil(p.b_c / p.dy);
        size_t ac_m  = (size_t) std::ceil(p.a_c / p.dx);
        size_t d_n  = (size_t) std::ceil(p.d / p.dy);
        size_t d_m  = (size_t) std::ceil(p.d / p.dx);
        size_t dc_n = (size_t) std::ceil(p.d_c / p.dy);
        size_t dc_m = (size_t) std::ceil(p.d_c / p.dx);

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
                bool is_cage_inner       = is_in_rect({ i, j }, { Y_n - yc_n - bc_n + dc_n, Y_n - yc_n + bc_n - dc_n, X_m + xc_m - ac_m + dc_m, X_m + xc_m + ac_m - dc_m });
                bool is_cage             = !is_cage_inner && is_in_rect({ i, j }, { Y_n - yc_n - bc_n, Y_n - yc_n + bc_n, X_m + xc_m - ac_m, X_m + xc_m + ac_m });
                bool is_cage_border      = is_cage && !is_in_rect({ i, j }, { Y_n - yc_n - bc_n + 1, Y_n - yc_n + bc_n - 1, X_m + xc_m - ac_m + 1, X_m + xc_m + ac_m - 1 })
                                        || !is_cage_inner && is_in_rect({ i, j }, { Y_n - yc_n - bc_n + dc_n - 1, Y_n - yc_n + bc_n - dc_n + 1, X_m + xc_m - ac_m + dc_m - 1, X_m + xc_m + ac_m - dc_m + 1 });

                if (is_border || is_capacitor_border || is_cage_border)
                    d.area_map[i][j]  = material::border;
                if (is_capacitor || is_cage || is_border)
                    d.area_map[i][j] |= material::metal;
                else
                    d.area_map[i][j]  = material::dielectr;

                if (is_capacitor && (j < Y_n)) d.fn[i][j] = p.q1;
                if (is_capacitor && (j > Y_n)) d.fn[i][j] = p.q2;
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

        d.coefs = [&] (size_t i, size_t j, chasing_dir dir,
                       const chasing_data & d,
                       const parameters & p,
                       std::vector < std::vector < double > > & T) -> chasing_coefs
        {
            // make ordinary A B C D coefficients
            auto make_normal_coefs = [&] (
                size_t i, size_t j, chasing_dir s,
                double eps) -> chasing_coefs
            {
                bool is_ext_i = (get_material_at(d, { (int) i - 1, ( int) j })
                                 | get_material_at(d, { (int) i + 1, ( int) j })) & material::ext;
                bool is_ext_j = (get_material_at(d, { (int) i, ( int) j - 1 })
                                 | get_material_at(d, { (int) i, ( int) j + 1 })) & material::ext;

                if (s == chasing_dir::i)
                {
                    return
                    {
                        - p.dt / 2. / p.dy / p.dy,
                        - p.dt / 2. / p.dy / p.dy,
                        1 + p.dt / 1. / p.dy / p.dy,
                        T[i][j] + p.dt / 2. *
                        (
                            (is_ext_j ? 0 : ((T[i][j + 1] + T[i][j - 1] - 2 * T[i][j]) / p.dx / p.dx))
                        ) + p.dt / 2. * d.fn[i][j] / eps
                    };
                }
                return
                {
                    - p.dt / 2. / p.dx / p.dx,
                    - p.dt / 2. / p.dx / p.dx,
                    1 + p.dt / 1. / p.dx / p.dx,
                    T[i][j] + p.dt / 2. *
                    (
                        (is_ext_i ? 0 : ((T[i + 1][j] + T[i - 1][j] - 2 * T[i][j]) / p.dy / p.dy))
                    ) + p.dt / 2. * d.fn[i][j] / eps
                };
            };

            material_t m = d.area_map[i][j];

            // deal with boundary conditions
            if (m & material::border)
            {
                if (m & material::metal)
                {
                    if (m & material::cap1)
                    {
                        return { 0, 0, 1, p.q1 };
                    }
                    else if (m & material::cap2)
                    {
                        return { 0, 0, 1, p.q2 };
                    }
                    else
                    {
                        return { 0, 0, 1, 0 };
                    }
                }
            }
            // deal with ordinary space
            else if (!(m & material::metal))
            {
                return make_normal_coefs(i, j, dir, p.eps);
            }

            return { 0, 0, 1, 0 };
        };
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

    inline static plot::painter_t make_system_painter(const parameters & params,
                                                      const chasing_data & d,
                                                      const std::vector < std::vector < double > > & T,
                                                      double & max_heat_value,
                                                      bool & draw_heat_map)
    {
        return [&] (CDC & dc, const plot::viewport & vp)
        {
            if (d.area_map.empty()) return;

            auto metal_brush  = plot::palette::brush(RGB(100, 100, 100));
            auto border_brush = plot::palette::brush(RGB( 50,  50,  50));

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
                            ((double)j - params.a * params.k - 0.5) * params.dx,
                            ((double)j - params.a * params.k + 0.5) * params.dx,
                            ((double)i - params.b * params.k - 0.5) * params.dy,
                            ((double)i - params.b * params.k + 0.5) * params.dy
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
