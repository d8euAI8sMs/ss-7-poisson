#pragma once

#include <vector>

#include <util/common/ptr.h>
#include <util/common/iterable.h>
#include <util/common/plot/palette.h>
#include <util/common/plot/drawable.h>
#include <util/common/plot/data_source.h>
#include <util/common/plot/custom_drawable.h>
#include <util/common/plot/utility.h>

namespace plot
{

    template
    <
        typename _container_t,
        typename _mapped_t = typename point < double >
    >
    class multilist_drawable : public drawable
    {

    public:

        using ptr_t = util::ptr_t < multilist_drawable > ;

        template < class ... T > static ptr_t create(T && ... t) { return util::create < typename ptr_t::element_type > (std::forward < T > (t) ...); }

    protected:

        using _iterator_t = typename _container_t::const_iterator;

    public:

        data_source_t < std::vector < _container_t > > data_source;
        util::data_mapper_t < _container_t, typename _container_t::const_iterator, _mapped_t > data_mapper;
        drawable::ptr_t  point_painter;
        palette::pen_ptr line_pen;

    public:

        multilist_drawable()
        {
        }

        multilist_drawable
        (
            data_source_t < std::vector < _container_t > > data_source,
            drawable::ptr_t  point_painter = custom_drawable::create(circle_painter()),
            palette::pen_ptr line_pen = palette::pen()
        )
            : data_source(std::move(data_source))
            , point_painter(std::move(point_painter))
            , line_pen(line_pen)
        {
        }

        virtual ~multilist_drawable()
        {
        }

    protected:

        virtual void do_draw(CDC &dc, const viewport &bounds) override
        {
            if ((!data_source) || (!data_mapper) && (!util::is_identity_mappable < _iterator_t, _mapped_t > ())) return;
            auto multidata = data_source(bounds);
            if (multidata->empty()) return;

            point < double > p;

            for each (auto & data in *multidata)
            {
                if (data.empty()) continue;

                dc.MoveTo(bounds.world_to_screen().xy(util::get_mapped_value_or_default(data_mapper, data, (_iterator_t)data.begin(), 0)));

                size_t idx = 0;
                _iterator_t end = data.end();

                for (_iterator_t it = data.begin(); it != end; ++it, ++idx)
                {
                    p = util::get_mapped_value_or_default(data_mapper, data, it, idx);
                    if (line_pen)
                    {
                        auto _old = dc.SelectObject(*line_pen);
                        dc.LineTo(bounds.world_to_screen().xy(p));
                        dc.SelectObject(_old);
                    }
                    if (point_painter)
                    {
                        point_painter->draw_at(dc, bounds, p);
                    }
                }
            }
        }
    };
}
