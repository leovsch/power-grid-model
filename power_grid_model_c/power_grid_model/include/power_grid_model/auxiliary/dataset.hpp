// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

// handle dataset and buffer related stuff

#include "../common/common.hpp"
#include "../common/counting_iterator.hpp"
#include "../common/exception.hpp"
#include "dataset_fwd.hpp"
#include "meta_data.hpp"

#include <boost/range.hpp>

#include <span>
#include <string_view>

namespace power_grid_model {

namespace meta_data {

template <dataset_type_tag T>
constexpr bool is_data_mutable_v = std::same_as<T, mutable_dataset_t> || std::same_as<T, writable_dataset_t>;
template <dataset_type_tag T> constexpr bool is_indptr_mutable_v = std::same_as<T, writable_dataset_t>;

static_assert(dataset_type_tag<const_dataset_t>);
static_assert(dataset_type_tag<mutable_dataset_t>);
static_assert(dataset_type_tag<writable_dataset_t>);
static_assert(!is_data_mutable_v<const_dataset_t>);
static_assert(is_data_mutable_v<mutable_dataset_t>);
static_assert(is_data_mutable_v<writable_dataset_t>);
static_assert(!is_indptr_mutable_v<const_dataset_t>);
static_assert(!is_indptr_mutable_v<mutable_dataset_t>);
static_assert(is_indptr_mutable_v<writable_dataset_t>);

struct ComponentInfo {
    MetaComponent const* component;
    // for non-uniform component, this is -1, we use indptr to describe the elements per scenario
    Idx elements_per_scenario;
    Idx total_elements;
};

struct DatasetInfo {
    bool is_batch;
    Idx batch_size; // for single dataset, the batch size is one
    MetaDataset const* dataset;
    std::vector<ComponentInfo> component_info;
};

template <typename T, dataset_type_tag dataset_type> class ColumnarAttributeRange {
  public:
    using Data = std::conditional_t<is_data_mutable_v<dataset_type>, void, void const>;

    class iterator;

    class Proxy {
      public:
        using value_type = T;

        Proxy() = default;
        Proxy(Idx idx, std::span<Data* const> data,
              std::span<std::reference_wrapper<MetaAttribute const> const> meta_attributes)
            : idx_{idx}, data_{data}, meta_attributes_{meta_attributes} {
            assert(data.size() == meta_attributes.size());
        }

        Proxy& operator=(value_type const& value)
            requires is_data_mutable_v<dataset_type>
        {
            for (Idx attribute_idx = 0; attribute_idx < meta_attributes_.size(); ++attribute_idx) {
                auto const& meta_attribute = get_meta_attribute(attribute_idx);
                char* buffer_ptr = reinterpret_cast<char*>(data_[attribute_idx]) + meta_attribute.size * idx_;
                meta_attribute.get_value(&value, buffer_ptr, 0);
            }
            return *this;
        }
        operator value_type() const { return get(); }
        value_type get() const {
            std::remove_const_t<value_type> result{};
            for (Idx attribute_idx = 0; attribute_idx < meta_attributes_.size(); ++attribute_idx) {
                auto const& meta_attribute = get_meta_attribute(attribute_idx);
                char const* buffer_ptr =
                    reinterpret_cast<char const*>(data_[attribute_idx]) + meta_attribute.size * idx_;
                meta_attribute.set_value(&result, buffer_ptr, 0);
            }
            return result;
        }

      private:
        friend class ColumnarAttributeRange;
        friend class iterator;

        MetaAttribute const& get_meta_attribute(Idx attribute_idx) const {
            return meta_attributes_[attribute_idx].get();
        }

        Idx idx_{};
        std::span<Data* const> data_{};
        std::span<std::reference_wrapper<MetaAttribute const> const> meta_attributes_;
    };

    class iterator : public boost::iterator_facade<iterator, T, boost::random_access_traversal_tag, Proxy, Idx> {
      public:
        using value_type = Proxy;

        iterator() = default;
        iterator(Idx idx, std::span<Data* const> data,
                 std::span<std::reference_wrapper<MetaAttribute const> const> meta_attributes)
            : current_{idx, data, meta_attributes} {}

      private:
        friend class boost::iterator_core_access;

        constexpr auto dereference() const -> value_type const { return *&current_; }
        constexpr auto dereference() -> value_type { return *&current_; }
        constexpr auto equal(iterator const& other) const {
            return current_.idx_ == other.current_.idx_ &&
                   current_.meta_attributes_.size() == other.current_.meta_attributes_.size() &&
                   current_.meta_attributes_.begin() == other.current_.meta_attributes_.begin();
        }
        constexpr auto distance_to(iterator const& other) const { return other.current_.idx_ - current_.idx_; }
        constexpr void increment() { ++current_.idx_; }
        constexpr void decrement() { --current_.idx_; }
        constexpr void advance(Idx n) { current_.idx_ += n; }

        Proxy current_;
    };

    ColumnarAttributeRange() = default;
    ColumnarAttributeRange(Idx size, std::vector<Data*> data,
                           std::span<std::reference_wrapper<MetaAttribute const> const> meta_attributes)
        : size_{size}, data_{std::move(data)}, meta_attributes_{meta_attributes} {
        assert(data_.size() == meta_attributes_.size());
    }
    ColumnarAttributeRange(ColumnarAttributeRange::iterator begin, ColumnarAttributeRange::iterator end)
        : size_{std::distance(begin, end)},
          data_{begin->data_.begin(), begin->data_.end()},
          meta_attributes_{begin->meta_attributes_} {
        assert(data_.size() == meta_attributes_.size());
        assert(begin + std::distance(begin, end) == end);
    }

    constexpr Idx size() const { return size_; }
    constexpr bool empty() const { return size_ == 0; }
    iterator begin() const { return get(0); }
    iterator end() const { return get(size_); }
    auto iter() const { return boost::make_iterator_range(begin(), end()); }
    auto operator[](Idx idx) const { return *get(idx); }

  private:
    iterator get(Idx idx) const { return iterator{idx, data_, meta_attributes_}; }

    Idx size_{};
    std::vector<Data*> data_{};
    std::span<std::reference_wrapper<MetaAttribute const> const> meta_attributes_;
};

template <typename T> using const_range_object = ColumnarAttributeRange<T, const_dataset_t>;
template <typename T> using mutable_range_object = ColumnarAttributeRange<T, mutable_dataset_t>;

template <dataset_type_tag dataset_type_> class Dataset {
    struct immutable_t {};
    struct mutable_t {};

  public:
    using dataset_type = dataset_type_;
    using Data = std::conditional_t<is_data_mutable_v<dataset_type>, void, void const>;
    using Indptr = std::conditional_t<is_indptr_mutable_v<dataset_type>, Idx, Idx const>;

    template <class StructType>
    using DataStruct = std::conditional_t<is_data_mutable_v<dataset_type>, StructType, StructType const>;

    struct AttributeBuffers {
        std::vector<Data*> data;
        std::vector<std::reference_wrapper<MetaAttribute const>> meta_attributes;
    };

    // for columnar buffers, Data* data is empty and attributes is filled
    // for uniform buffers, indptr is empty
    struct Buffer {
        Data* data;
        AttributeBuffers attributes;
        std::span<Indptr> indptr;
    };

    template <class StructType>
    using RangeObject = std::conditional_t<is_data_mutable_v<dataset_type>, mutable_range_object<StructType>,
                                           const_range_object<StructType>>;

    static constexpr Idx invalid_index{-1};

    Dataset(bool is_batch, Idx batch_size, std::string_view dataset_name, MetaData const& meta_data)
        : meta_data_{&meta_data},
          dataset_info_{.is_batch = is_batch,
                        .batch_size = batch_size,
                        .dataset = &meta_data.get_dataset(dataset_name),
                        .component_info = {}} {
        if (dataset_info_.batch_size < 0) {
            throw DatasetError{"Batch size cannot be negative!\n"};
        }
        if (!dataset_info_.is_batch && (dataset_info_.batch_size != 1)) {
            throw DatasetError{"For non-batch dataset, batch size should be one!\n"};
        }
    }

    // implicit conversion constructor from writable to mutable, from writable and mutable to const
    template <dataset_type_tag other_dataset_type>
        requires((is_data_mutable_v<other_dataset_type> && !is_data_mutable_v<dataset_type>) ||
                 (is_indptr_mutable_v<other_dataset_type> && !is_indptr_mutable_v<dataset_type>))
    Dataset(Dataset<other_dataset_type> const& other)
        : meta_data_{&other.meta_data()}, dataset_info_{other.get_description()} {
        for (Idx i{}; i != other.n_components(); ++i) {
            auto const& buffer = other.get_buffer(i);
            buffers_.push_back(Buffer{.data = buffer.data, .indptr = buffer.indptr});
        }
    }

    MetaData const& meta_data() const { return *meta_data_; }
    bool empty() const { return dataset_info_.component_info.empty(); }
    bool is_batch() const { return dataset_info_.is_batch; }
    Idx batch_size() const { return dataset_info_.batch_size; }
    MetaDataset const& dataset() const { return *dataset_info_.dataset; }
    Idx n_components() const { return static_cast<Idx>(buffers_.size()); }
    DatasetInfo const& get_description() const { return dataset_info_; }
    ComponentInfo const& get_component_info(Idx i) const { return dataset_info_.component_info[i]; }
    Buffer const& get_buffer(std::string_view component) const { return get_buffer(find_component(component, true)); }
    Buffer const& get_buffer(Idx i) const { return buffers_[i]; }

    Idx find_component(std::string_view component, bool required = false) const {
        auto const found = std::ranges::find_if(dataset_info_.component_info, [component](ComponentInfo const& x) {
            return x.component->name == component;
        });
        if (found == dataset_info_.component_info.cend()) {
            if (required) {
                using namespace std::string_literals;
                throw DatasetError{"Cannot find component '"s + std::string{component} + "'!\n"s};
            }
            return invalid_index;
        }
        return std::distance(dataset_info_.component_info.cbegin(), found);
    }
    bool contains_component(std::string_view component) const { return find_component(component) >= 0; }

    ComponentInfo const& get_component_info(std::string_view component) const {
        return dataset_info_.component_info[find_component(component, true)];
    }

    void add_component_info(std::string_view component, Idx elements_per_scenario, Idx total_elements)
        requires is_indptr_mutable_v<dataset_type>
    {
        add_component_info_impl(component, elements_per_scenario, total_elements);
    }

    void add_buffer(std::string_view component, std::integral auto elements_per_scenario_,
                    std::integral auto total_elements_, Indptr* indptr, Data* data)
        requires(!is_indptr_mutable_v<dataset_type>)
    {
        auto const elements_per_scenario = static_cast<Idx>(elements_per_scenario_);
        auto const total_elements = static_cast<Idx>(total_elements_);
        check_non_uniform_integrity<immutable_t>(elements_per_scenario, total_elements, indptr);
        add_component_info_impl(component, elements_per_scenario, total_elements);
        buffers_.back().data = data;
        if (indptr) {
            buffers_.back().indptr = get_indptr_span(indptr);
        } else {
            buffers_.back().indptr = {};
        }
    }

    void set_buffer(std::string_view component, Indptr* indptr, Data* data)
        requires is_indptr_mutable_v<dataset_type>
    {
        Idx const idx = find_component(component, true);
        ComponentInfo const& info = dataset_info_.component_info[idx];
        check_non_uniform_integrity<mutable_t>(info.elements_per_scenario, info.total_elements, indptr);
        buffers_[idx].data = data;
        if (indptr) {
            buffers_[idx].indptr = get_indptr_span(indptr);
        } else {
            buffers_[idx].indptr = {};
        }
    }

    void add_attribute_buffer(std::string_view component, std::string_view attribute, Data* data) {
        Idx const idx = find_component(component, true);
        Buffer& buffer = buffers_[idx];
        if (buffer.data != nullptr) {
            throw DatasetError{"Cannot add attribute buffers to row-based dataset!\n"};
        }
        if (std::ranges::find_if(buffer.attributes.meta_attributes, [&attribute](MetaAttribute const& meta_attribute) {
                return meta_attribute.name == attribute;
            }) != buffer.attributes.meta_attributes.end()) {
            throw DatasetError{"Cannot have duplicated attribute buffers!\n"};
        }
        buffer.attributes.data.emplace_back(data);
        buffer.attributes.meta_attributes.emplace_back(
            dataset_info_.component_info[idx].component->get_attribute(attribute));
    }

    // get buffer by component type
    template <template <class> class type_getter, class ComponentType,
              class StructType = DataStruct<typename type_getter<ComponentType>::type>>
    std::span<StructType> get_buffer_span(Idx scenario = invalid_index) const {
        assert(scenario < batch_size());

        if (!is_batch() && scenario > 0) {
            throw DatasetError{"Cannot export a single dataset with specified scenario\n"};
        }

        Idx const idx = find_component(ComponentType::name, false);
        return get_buffer_span_impl<StructType>(scenario, idx);
    }
    template <template <class> class type_getter, class ComponentType,
              class StructType = DataStruct<typename type_getter<ComponentType>::type>>
    RangeObject<StructType> get_columnar_buffer_span(Idx scenario = invalid_index) const {
        assert(scenario < batch_size());

        if (!is_batch() && scenario > 0) {
            throw DatasetError{"Cannot export a single dataset with specified scenario\n"};
        }

        Idx const idx = find_component(ComponentType::name, false);
        return get_columnar_buffer_span_impl<StructType>(scenario, idx);
    }

    // get buffer by component type for all scenarios in vector span
    template <template <class> class type_getter, class ComponentType,
              class StructType = DataStruct<typename type_getter<ComponentType>::type>>
    std::vector<std::span<StructType>> get_buffer_span_all_scenarios() const {
        Idx const idx = find_component(ComponentType::name, false);
        std::vector<std::span<StructType>> result(batch_size());
        for (Idx scenario{}; scenario != batch_size(); scenario++) {
            result[scenario] = get_buffer_span_impl<StructType>(scenario, idx);
        }
        return result;
    }
    template <template <class> class type_getter, class ComponentType,
              class StructType = DataStruct<typename type_getter<ComponentType>::type>>
    std::vector<RangeObject<StructType>> get_columnar_buffer_span_all_scenarios() const {
        Idx const idx = find_component(ComponentType::name, false);
        std::vector<RangeObject<StructType>> result(batch_size());
        for (Idx scenario{}; scenario != batch_size(); scenario++) {
            result[scenario] = get_columnar_buffer_span_impl<StructType>(scenario, idx);
        }
        return result;
    }

    // get individual dataset from batch
    Dataset get_individual_scenario(Idx scenario)
        requires(!is_indptr_mutable_v<dataset_type>)
    {
        assert(0 <= scenario && scenario < batch_size());

        Dataset result{false, 1, dataset().name, meta_data()};
        for (Idx i{}; i != n_components(); ++i) {
            auto const& buffer = get_buffer(i);
            auto const& component_info = get_component_info(i);
            Idx size = component_info.elements_per_scenario >= 0
                           ? component_info.elements_per_scenario
                           : buffer.indptr[scenario + 1] - buffer.indptr[scenario];
            Data* data = component_info.elements_per_scenario >= 0
                             ? component_info.component->advance_ptr(buffer.data, size * scenario)
                             : component_info.component->advance_ptr(buffer.data, buffer.indptr[scenario]);
            result.add_buffer(component_info.component->name, size, size, nullptr, data);
        }
        return result;
    }

  private:
    MetaData const* meta_data_;
    DatasetInfo dataset_info_;
    std::vector<Buffer> buffers_;

    std::span<Indptr> get_indptr_span(Indptr* indptr) const {
        return std::span{indptr, static_cast<size_t>(batch_size() + 1)};
    }

    void check_uniform_integrity(Idx elements_per_scenario, Idx total_elements) {
        if ((elements_per_scenario >= 0) && (elements_per_scenario * batch_size() != total_elements)) {
            throw DatasetError{
                "For a uniform buffer, total_elements should be equal to elements_per_scenario * batch_size!\n"};
        }
    }

    template <typename check_indptr_content>
        requires std::same_as<check_indptr_content, mutable_t> || std::same_as<check_indptr_content, immutable_t>
    void check_non_uniform_integrity(Idx elements_per_scenario, Idx total_elements, Indptr* indptr) {
        if (elements_per_scenario < 0) {
            if (!indptr) {
                throw DatasetError{"For a non-uniform buffer, indptr should be supplied!\n"};
            }
            if constexpr (std::same_as<check_indptr_content, immutable_t>) {
                if (indptr[0] != 0 || indptr[batch_size()] != total_elements) {
                    throw DatasetError{
                        "For a non-uniform buffer, indptr should begin with 0 and end with total_elements!\n"};
                }
            }
        } else if (indptr) {
            throw DatasetError{"For a uniform buffer, indptr should be nullptr!\n"};
        }
    }

    void add_component_info_impl(std::string_view component, Idx elements_per_scenario, Idx total_elements) {
        if (find_component(component) >= 0) {
            throw DatasetError{"Cannot have duplicated components!\n"};
        }
        check_uniform_integrity(elements_per_scenario, total_elements);
        dataset_info_.component_info.push_back(
            {&dataset_info_.dataset->get_component(component), elements_per_scenario, total_elements});
        buffers_.push_back(Buffer{});
    }

    template <class StructType> std::span<StructType> get_buffer_span_impl(Idx scenario, Idx component_idx) const {
        // return empty span if the component does not exist
        if (component_idx < 0) {
            return {};
        }
        // return span based on uniform or non-uniform buffer
        ComponentInfo const& info = dataset_info_.component_info[component_idx];
        Buffer const& buffer = buffers_[component_idx];
        auto const ptr = reinterpret_cast<StructType*>(buffer.data);
        if (scenario < 0) {
            return std::span<StructType>{ptr, ptr + info.total_elements};
        }
        if (info.elements_per_scenario < 0) {
            return std::span<StructType>{ptr + buffer.indptr[scenario], ptr + buffer.indptr[scenario + 1]};
        }
        return std::span<StructType>{ptr + info.elements_per_scenario * scenario,
                                     ptr + info.elements_per_scenario * (scenario + 1)};
    }

    template <class StructType>
    RangeObject<StructType> get_columnar_buffer_span_impl(Idx scenario, Idx component_idx) const {
        // return empty span if the component does not exist
        if (component_idx < 0) {
            return {};
        }
        // return span based on uniform or non-uniform buffer
        ComponentInfo const& info = dataset_info_.component_info[component_idx];
        Buffer const& buffer = buffers_[component_idx];
        assert(buffer.data == nullptr);

        RangeObject<StructType> total_range{info.total_elements, buffer.attributes.data,
                                            buffer.attributes.meta_attributes};
        if (scenario < 0) {
            return total_range;
        }
        // if (info.elements_per_scenario < 0) {
        //     return RangeObject<StructType>{ptr + buffer.indptr[scenario], ptr + buffer.indptr[scenario + 1]};
        // }
        return RangeObject<StructType>{total_range.begin() + scenario * info.elements_per_scenario,
                                       total_range.begin() + info.elements_per_scenario * (scenario + 1)};
    }
};

} // namespace meta_data

template <dataset_type_tag dataset_type> using Dataset = meta_data::Dataset<dataset_type>;
using ConstDataset = Dataset<const_dataset_t>;
using MutableDataset = Dataset<mutable_dataset_t>;
using WritableDataset = Dataset<writable_dataset_t>;

} // namespace power_grid_model
