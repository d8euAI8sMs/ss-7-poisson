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

        triangle () = default;

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

    bool symmetric = false;

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

    double range_of_T;

    if ((global_max_T > 0) && (global_min_T < 0)
        && (std::abs(global_max_T + global_min_T) / global_max_T) < 0.1)
    {
        symmetric = true;
        global_min_T = 0;
    }

    range_of_T = min(global_max_T - global_min_T, + dT * max_isolines);

    for (size_t i = 0; i < n - 1; ++i)
    {
        visited[i].resize(m);
        for (size_t j = 0; j < m - 1; ++j)
        {
            if (!stencil(i, j)) continue;

            double max_T = max(max(T[i][j], T[i][j + 1]), max(T[i + 1][j], T[i + 1][j + 1]));
            double min_T = min(min(T[i][j], T[i][j + 1]), min(T[i + 1][j], T[i + 1][j + 1]));

            if (!symmetric) max_T = min(global_min_T + range_of_T, max_T);
            else {
                max_T = min(range_of_T, max_T);
                min_T = max(-range_of_T, min_T);
            }
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

inline static plot::point < double > gradiend_at
(
    const std::vector < std::vector < double > > & T,
    size_t n, size_t m,
    const triangle & t
)
{
    size_t i = t.i, j = t.j;

    double ni, nj;

    if (t.o == 0)
    {
        ni = T[i][j + 1] - T[i][j];
        nj = T[i + 1][j] - T[i][j];
    }
    else
    {
        ni = - (T[i + 1][j] - T[i + 1][j + 1]);
        nj = - (T[i][j + 1] - T[i + 1][j + 1]);
    }

    double norm = std::sqrt(ni * ni + nj * nj);

    if (!isfinite(ni / norm) || !isfinite(nj / norm))
    {
        return { 1, 0 };
    }

    return { ni / norm, nj / norm };
}

inline static triangle get_containing_triangle(plot::point < double > p)
{
    size_t i = (size_t) std::trunc(p.y);
    size_t j = (size_t) std::trunc(p.x);

    return { i, j, ((p.y - (double)j) < (1 - (p.x - (double) i))) ? 0 : 1u };
}

inline static bool next_point
(
    const std::vector < std::vector < double > > & T,
    const std::vector < std::vector < std::array < owning_type, 2 > > > & visited,
    size_t n, size_t m,
    std::pair < triangle, plot::point < double > > & cur,
    bool inv
)
{
    auto & t = std::get<0>(cur);
    auto & p = std::get<1>(cur);

    auto ngrad = gradiend_at(T, n, m, t);

    const double d = (inv ? -0.1 : 0.1);

    auto c = t;

    do
    {
        p = { p.x + d * ngrad.x, p.y + d * ngrad.y };
        t = get_containing_triangle(p);
    } while ((c.i == t.i) && (c.j == t.j) && (c.o == t.o));

    if (visited[t.i][t.j][t.o] == owning_type::our) return false;

    return true;
}

void model::find_field_lines
(
    const std::vector < std::vector < double > > & T,
    std::vector < std::vector < plot::point < double > > > & out,
    size_t n, size_t m,
    const parameters & p,
    stencil_fn stencil,
    std::vector < plot::point < double > > hint
)
{
    std::vector < std::vector < std::array < owning_type, 2 > > > visited(
        n, std::vector < std::array < owning_type, 2 > >(m));

    std::vector < triangle > path;
    std::vector < plot::point < double > > isoline;

    for each (auto ph in hint)
    {
        plot::point < size_t > p0 = { (size_t) std::trunc(ph.y), (size_t) std::trunc(ph.x) };

        if
        (
            (p0.x <= 0)
         || (p0.y <= 0)
         || (p0.x >= (n - 1))
         || (p0.y >= (m - 1))
         || !stencil(p0.x, p0.y)
        ) continue;

        std::pair < triangle, plot::point < double > > cur;

        auto & t = std::get<0>(cur);

        bool success;

        bool inv = true;

        // two iterations (in forward and backward direction)

        do
        {
            inv = !inv;

            cur =
            {
                get_containing_triangle(ph),
                ph
            };

            visited[t.i][t.j][t.o] = owning_type::our;

            path.push_back(std::get<0>(cur));
            isoline.emplace_back(
                - p.a * p.k + std::get<1>(cur).x * p.dx,
                - p.b * p.k + std::get<1>(cur).y * p.dy);

            do
            {
                success = next_point(T, visited, n, m, cur, inv);
                if (success)
                {
                    visited[t.i][t.j][t.o] = owning_type::our;
                    path.push_back(std::get<0>(cur));
                    isoline.emplace_back(
                        - p.a * p.k + std::get<1>(cur).x * p.dx,
                        - p.b * p.k + std::get<1>(cur).y * p.dy);
                }
            }
            while ((t.i > 0) && (t.j > 0) && (t.i < (n - 1)) && (t.j < (m - 1))
                     && success && stencil(t.i, t.j));

            out.push_back(std::move(isoline));

            for each (auto & t in path)
            {
                visited[t.i][t.j][t.o] = owning_type::their;
            }

            path.clear();
        }
        while (!inv);
    }
}
