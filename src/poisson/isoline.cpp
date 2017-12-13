#include "stdafx.h"

#include <vector>
#include <array>
#include <set>

#include "model.h"

using namespace model;

namespace
{

    /**
     * The `triangle` struct describes one of
     * the following triangle objects:
     *
     *       (i,j)   (i,j+1)
     *           +---+
     *           |  /|
     * o == 0 -> | / | <- o == 1
     *           |/  |
     *           +---+
     *     (i+1,j)   (i+1,j+1)
     */

    struct triangle
    {
        size_t i, j, o;

        triangle (size_t i, size_t j, size_t o)
            : i(i), j(j), o(o)
        {
        }
    };

    enum class side : size_t
    {
        left, top, right, bottom, diagonal
    };

    enum class owning_type : size_t
    {
        no, our, their
    };

    inline bool is_in_range(double c, double a, double b)
    {
        return (((a <= c) && (c < b)) || ((b <= c) && (c < a)));
    }

    inline plot::point < double > get_pos_in_range
    (
        const parameters & p,
        double T0,
        const std::vector < std::vector < double > > & T,
        plot::point < size_t > a,
        plot::point < size_t > b
    )
    {
        if (a.x == b.x)
        {
            return
            {
                - p.a * p.k + ((T0 - T[a.x][a.y]) / (T[b.x][b.y] - T[a.x][a.y]) * ((int) b.y - (int) a.y) + a.y) * p.dx,
                - p.b * p.k + a.x * p.dy
            }; // { x, y }
        }
        else if (a.y == b.y)
        {
            return
            {
                - p.a * p.k + a.y * p.dx,
                - p.b * p.k + ((T0 - T[a.x][a.y]) / (T[b.x][b.y] - T[a.x][a.y]) * ((int) b.x - (int) a.x) + a.x) * p.dy
            }; // { x, y }
        }
        else
        {
            return
            {
                - p.a * p.k + ((T0 - T[a.x][a.y]) / (T[b.x][b.y] - T[a.x][a.y]) * ((int) b.y - (int) a.y) + a.y) * p.dx,
                - p.b * p.k + ((T0 - T[a.x][a.y]) / (T[b.x][b.y] - T[a.x][a.y]) * ((int) b.x - (int) a.x) + a.x) * p.dy
            }; // { x, y }
        }
    }

    /* with respect to t1 */
    inline side get_common_side
    (
        const triangle & t1,
        const triangle & t2
    )
    {
        if ((t1.i == t2.i) && (t1.j == t2.j))
        {
            return side::diagonal;
        }
        else if (t1.j == t2.j + 1)
        {
            return side::left;
        }
        else if (t1.j + 1 == t2.j)
        {
            return side::right;
        }
        else if (t1.i == t2.i + 1)
        {
            return side::top;
        }
        else if (t1.i + 1 == t2.i)
        {
            return side::bottom;
        }

        throw std::exception("unexpected");
    }

    /* get side intersected by isoline except the given side */
    inline side get_side_with_intersection
    (
        const triangle & t,
        side except,
        double T0,
        const std::vector < std::vector < double > > & T
    )
    {
        size_t i = t.i, j = t.j;
        bool has_left   = is_in_range(T0, T[i][j], T[i + 1][j]);
        bool has_top    = is_in_range(T0, T[i][j], T[i][j + 1]);
        bool has_right  = is_in_range(T0, T[i][j + 1], T[i + 1][j + 1]);
        bool has_bottom = is_in_range(T0, T[i + 1][j], T[i + 1][j + 1]);
        bool has_diag   = is_in_range(T0, T[i + 1][j], T[i][j + 1]);

        if (t.o == 0)
        {
            if (has_left && (except != side::left)) return side::left;
            if (has_top && (except != side::top)) return side::top;
        }
        else
        {
            if (has_right && (except != side::right)) return side::right;
            if (has_bottom && (except != side::bottom)) return side::bottom;
        }
        if (has_diag && (except != side::diagonal)) return side::diagonal;

        return except;
    }

    inline void side_to_point
    (
        const triangle & t,
        side s,
        plot::point < size_t > & p1,
        plot::point < size_t > & p2
    )
    {
        switch (s)
        {
        case side::left:
            p1 = { t.i, t.j };
            p2 = { t.i + 1, t.j };
            break;
        case side::top:
            p1 = { t.i, t.j };
            p2 = { t.i, t.j + 1 };
            break;
        case side::right:
            p1 = { t.i, t.j + 1 };
            p2 = { t.i + 1, t.j + 1 };
            break;
        case side::bottom:
            p1 = { t.i + 1, t.j };
            p2 = { t.i + 1, t.j + 1 };
            break;
        case side::diagonal:
            p1 = { t.i + 1, t.j };
            p2 = { t.i, t.j + 1 };
            break;
        default:
            throw std::exception("unexpected");
        }
    }
}

void model::find_isolines
(
    const std::vector < std::vector < double > > & T,
    double dT,
    std::vector < std::vector < plot::point < double > > > & out,
    size_t n, size_t m,
    const parameters & p,
    stencil_fn stencil,
    size_t max_isolines,
    size_t max_points_in_stack
)
{
    /* T0_k index => triangle */
    std::vector < std::pair < int, triangle > > stack;

    /* [i][j][o] => < owning type, owner > */
    std::vector
    <
        std::vector
        <
            std::array < std::pair < owning_type, std::set < int > >, 2 >
        >
    >
    visited(n);

    std::vector < triangle > path;

    /* Search for triangles intersected by isolines */

    double global_max_T = (std::numeric_limits < double > :: lowest)();
    double global_min_T = (std::numeric_limits < double > :: max)();

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            if (!stencil(i, j)) continue;
            if (global_max_T < T[i][j]) global_max_T = T[i][j];
            if (global_min_T > T[i][j]) global_min_T = T[i][j];
        }
    }

    global_max_T = min(global_max_T, global_min_T + dT * max_isolines);

    for (size_t i = 0; i < n - 1; ++i)
    {
        visited[i].resize(m);
        for (size_t j = 0; j < m - 1; ++j)
        {
            if (!stencil(i, j)) continue;

            double max_T = min(global_max_T, max(max(T[i][j], T[i][j + 1]), max(T[i + 1][j], T[i + 1][j + 1])));
            double min_T = min(min(T[i][j], T[i][j + 1]), min(T[i + 1][j], T[i + 1][j + 1]));
            if (min_T < max_T)
            {
                int a = std::ceil(max_T / dT);
                int b = std::floor(min_T / dT);
                size_t r = a - b;
                for (size_t k = 0; k < r; ++k)
                {
                    double T0 = (b + (int) k) * dT;

                    /* o == 0 type */
                    if (stencil(i + 1, j) && is_in_range(T0, T[i][j], T[i + 1][j]) ||
                        stencil(i, j + 1) && is_in_range(T0, T[i][j], T[i][j + 1]))
                    {
                        stack.emplace_back((b + (int) k), triangle { i, j, 0 });
                    }

                    /* o == 1 type */
                    if (stencil(i + 1, j + 1) && stencil(i, j + 1) && is_in_range(T0, T[i][j + 1], T[i + 1][j + 1]) ||
                        stencil(i + 1, j) && stencil(i + 1, j + 1) && is_in_range(T0, T[i + 1][j], T[i + 1][j + 1]))
                    {
                        stack.emplace_back((b + (int) k), triangle { i, j, 1 });
                    }

                    if (stack.size() > max_points_in_stack) return;
                }
            }
        }
    }

    /* Trace an isoline at the top of stack */

    size_t current_contour = 0;

    while (!stack.empty())
    {
        size_t   Tk = std::get < 0 > (stack.back());
        double   T0 = std::get < 0 > (stack.back()) * dT;
        triangle t  = std::get < 1 > (stack.back());
        size_t  &i = t.i,
                &j = t.j;

        stack.pop_back();

        if (std::get < 0 > (visited[i][j][t.o]) != owning_type::no)
        {
            if (std::get < 1 > (visited[i][j][t.o]).find(Tk)
                != std::get < 1 > (visited[i][j][t.o]).end())
            {
                continue;
            }
        }

        for (size_t c = current_contour; c == current_contour;)
        {

            std::get < 0 > (visited[i][j][t.o]) = owning_type::our;
            std::get < 1 > (visited[i][j][t.o]).insert(Tk);

            path.push_back(t);

            if (t.o == 0)
            {
                bool has_left = stencil(i + 1, j) && is_in_range(T0, T[i][j], T[i + 1][j]);
                bool has_top  = stencil(i, j + 1) && is_in_range(T0, T[i][j], T[i][j + 1]);
                bool has_diag = stencil(i + 1, j) && stencil(i, j + 1) && is_in_range(T0, T[i + 1][j], T[i][j + 1]);

                if (stencil(i, (int) j - 1) && has_left && (std::get < 0 > (visited[i][j - 1][1]) != owning_type::our))
                {
                    t = { i, j - 1, 1 };
                }
                else if (stencil((int) i - 1, j) && has_top && (std::get < 0 > (visited[i - 1][j][1]) != owning_type::our))
                {
                    t = { i - 1, j, 1 };
                }
                else if (has_diag && (std::get < 0 > (visited[i][j][1]) != owning_type::our))
                {
                    t = { i, j, 1 };
                }
                else
                {
                    ++current_contour;
                }
            }
            else
            {
                bool has_right  = stencil(i, j + 1) && stencil(i + 1, j + 1) && is_in_range(T0, T[i][j + 1], T[i + 1][j + 1]);
                bool has_bottom = stencil(i + 1, j) && stencil(i + 1, j + 1) && is_in_range(T0, T[i + 1][j], T[i + 1][j + 1]);
                bool has_diag   = stencil(i + 1, j) && stencil(i, j + 1) && is_in_range(T0, T[i + 1][j], T[i][j + 1]);

                if (stencil(i, j + 2) && has_right && (std::get < 0 > (visited[i][j + 1][0]) != owning_type::our))
                {
                    t = { i, j + 1, 0 };
                }
                else if (stencil(i + 2, j) && has_bottom && (std::get < 0 > (visited[i + 1][j][0]) != owning_type::our))
                {
                    t = { i + 1, j, 0 };
                }
                else if (has_diag && (std::get < 0 > (visited[i][j][0]) != owning_type::our))
                {
                    t = { i, j, 0 };
                }
                else
                {
                    ++current_contour;
                }
            }
        }

        if (!path.empty())
        {
            // todo this, set visited to they
            std::vector < plot::point < double > > isoline(path.size() + 1);

            triangle t0 = path[0];
            size_t  &i0 = t0.i,
                    &j0 = t0.j;

            std::get < 0 > (visited[i0][j0][t0.o]) = owning_type::their;

            plot::point < size_t > p1, p2;

            for (size_t k = 1; k < path.size(); ++k)
            {
                t = path[k];

                std::get < 0 > (visited[i][j][t.o]) = owning_type::their;

                if (k == 1)
                {
                    side common = get_common_side(t0, t);
                    side other  = get_side_with_intersection(t0, common, T0, T);
                    side_to_point(t0, other, p1, p2);
                    isoline[0] = get_pos_in_range(p, T0, T, p1, p2);
                    side_to_point(t0, common, p1, p2);
                    isoline[1] = get_pos_in_range(p, T0, T, p1, p2);
                }

                side common = get_common_side(t, t0);
                side other  = get_side_with_intersection(t, common, T0, T);
                side_to_point(t, other, p1, p2);
                isoline[k + 1] = get_pos_in_range(p, T0, T, p1, p2);

                t0 = t;
            }

            out.push_back(std::move(isoline));

            path.clear();
        }
    }
}
